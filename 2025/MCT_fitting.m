%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 코드 : 특정 MCT 번호만 입력받아,
%        2-RC 모델 (R0, R1, R2, tau1, tau2) 파라미터 (5개) 피팅 예시
%        - packVoltage -> cellVoltage (/192)
%        - packCurrent -> cellCurrent (/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 코드1에서 저장한 결과 불러오기
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

% Q_batt 후보 (예: 50 ~ 60 사이 10개 구간)
Q_batt_candidates = 55.6;%linspace(50,60,10);  % [Ah]라고 가정

%% 3) fmincon 옵션 및 파라미터 초기값/제약조건 설정
%  5개 파라미터: X = [R0, R1, R2, tau1, tau2]
%    - tau1 = R1*C1
%    - tau2 = R2*C2

x0 = [0.0038,  ... R0
      0.02,   ... R1
      0.02,   ... R2
      3,     ... tau1 
      25];  ... tau2 

lb = [0, 0, 0, 0, 0];
ub = [inf, inf, inf, 6, 180];

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                       'MaxIterations',1000,'MaxFunctionEvaluations',5000);

%% 4) 해당 MCT 데이터
dataMCT = mctCellData{mctNumber};

time_s      = dataMCT.Time_s;
packVoltage = dataMCT.PackVoltage_V;
packCurrent = dataMCT.Current_A;   % +가 방전

% (a) "셀 전압"으로 변환 (192 직렬)
cellVoltage_meas = packVoltage / 192;

% (b) "셀 전류"로 변환 (2 병렬)
cellCurrent = packCurrent / 2;

%% (중요) 초기 SOC 계산 (OCV 테이블 기반 보간)
cellVoltage_init = mean(cellVoltage_meas(1:5));  
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear', 'extrap');

% 시계열 간격
dt = [0; diff(time_s)];

%% 5) Q_batt 후보별로 (R0, R1, R2, tau1, tau2) 최적화
rmseArray = zeros(length(Q_batt_candidates), 1);
paramArray = zeros(length(Q_batt_candidates), 5);  % [R0, R1, R2, tau1, tau2]

for iQ = 1:length(Q_batt_candidates)
    Qbatt = Q_batt_candidates(iQ);  % [Ah]
    
    % (A) SOC(t) 계산
    %     cellCurrent가 +면 방전
    charge_integral = cumtrapz(time_s, cellCurrent); % [A·s]
    SOC_t = SOC0 - (charge_integral / (Qbatt * 3600)) * 100;
    
    % (B) fmincon으로 최적화 (RMSE 최소화)
    costFunc = @(X) computeRMSE_2RC_tau(X, SOC_t, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
    [xOpt, fVal] = fmincon(costFunc, x0, [], [], [], [], lb, ub, [], options);
    
    rmseArray(iQ) = fVal;
    paramArray(iQ, :) = xOpt;  % (R0, R1, R2, tau1, tau2)
end

% (C) 최소 RMSE 지점
[bestRMSE, idxMin] = min(rmseArray);
bestQbatt = Q_batt_candidates(idxMin);
bestParams = paramArray(idxMin, :);  % [R0, R1, R2, tau1, tau2]
bestR0   = bestParams(1);
bestR1   = bestParams(2);
bestR2   = bestParams(3);
bestTau1 = bestParams(4);
bestTau2 = bestParams(5);

%% 6) 결과 출력
fprintf('\n=== MCT-%d 결과 (2-RC, 셀 전압; R0,R1,R2,tau1,tau2) ===\n', mctNumber);
fprintf('Q_batt 후보: '); fprintf('%.1f ', Q_batt_candidates); fprintf('\n');
fprintf('RMSE들   : '); fprintf('%.4f ', rmseArray); fprintf('\n');
fprintf('--> 최소 RMSE = %.4f @ Q_batt = %.2f [Ah]\n', bestRMSE, bestQbatt);
fprintf('    (R0 = %.4f, R1 = %.4f, R2 = %.4f, tau1 = %.1f, tau2 = %.1f)\n', ...
        bestR0, bestR1, bestR2, bestTau1, bestTau2);

%% 6.1) FittingResult Table 생성 (한눈에 보기 쉽게)
FittingResult = table(Q_batt_candidates(:), rmseArray, ...
                     paramArray(:, 1), paramArray(:, 2), paramArray(:, 3), ...
                     paramArray(:, 4), paramArray(:, 5), ...
                     'VariableNames', {'Q_batt', 'RMSE', 'R0', 'R1', 'R2', 'tau1', 'tau2'});
disp('--- FittingResult Table ---');
disp(FittingResult);

%% 7) Q_batt vs RMSE 플롯
figure('Name', ['MCT-' num2str(mctNumber) ' : RMSE vs Q_batt (2-RC)'], 'NumberTitle', 'off');
plot(Q_batt_candidates, rmseArray, 'o--', 'LineWidth', 1.2);
hold on;
plot(bestQbatt, bestRMSE, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Q_{batt} [Ah]');
ylabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : 2-RC 모델 RMSE vs Capacity']);
grid on;

%% 8) 최적 파라미터로 전압 비교 플롯
% SOC 재계산
charge_integral_best = cumtrapz(time_s, cellCurrent);
SOC_t_best = SOC0 - (charge_integral_best / (bestQbatt * 3600)) * 100;

% 모델 전압 계산
V_estBest = modelVoltage_2RC_tau(bestParams, SOC_t_best, cellCurrent, dt, socOCV, ocvCellVoltage);

figure('Name', ['MCT-' num2str(mctNumber) ' : Voltage Comparison (2-RC)'], 'NumberTitle', 'off');
subplot(2,1,1);
plot(time_s, cellVoltage_meas, '-b', 'LineWidth', 1.2); hold on;
plot(time_s, V_estBest, '--r', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Cell Voltage (V)');
legend('Measured', 'Model', 'Location', 'best');
title(['MCT-' num2str(mctNumber) ...
       ' 전압 비교 (2-RC, Q_{batt} = ' num2str(bestQbatt, '%.2f') 'Ah)']);
grid on;

subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_estBest, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
title('Error = (Measured - Model)');
grid on;

function cost = computeRMSE_2RC_tau(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    % X: [R0, R1, R2, tau1, tau2]
    V_est = modelVoltage_2RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost = sqrt(mean((V_est - V_meas).^2));
end

function V_est = modelVoltage_2RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % X: [R0, R1, R2, tau1, tau2]
    R0   = X(1);
    R1   = X(2);
    R2   = X(3);
    tau1 = X(4);
    tau2 = X(5);
    
    % 2개의 RC(이제는 tau1, tau2) 전압 항 (초기값)
    Vrc1 = 0;
    Vrc2 = 0;
    
    N = length(SOC_t);
    V_est = zeros(N, 1);
    
    for k = 1:N
        % (1) OCV 추출 (SOC에 따른 OCV 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');
        
        % (2) R0에 의한 전압 강하
        IR0 = R0 * I_cell(k);
        
        % (3) RC1, RC2 요소 업데이트
        if k > 1
            % alpha1 = exp(-dt(k)/(R1*C1)) => exp(-dt(k)/tau1)
            alpha1 = exp(-dt(k) / tau1);
            alpha2 = exp(-dt(k) / tau2);
            
            Vrc1 = Vrc1 * alpha1 + R1 * (1 - alpha1) * I_cell(k);
            Vrc2 = Vrc2 * alpha2 + R2 * (1 - alpha2) * I_cell(k);
        end
        
        % (4) 최종 모델 전압
        V_est(k) = OCV_now - IR0 - Vrc1 - Vrc2;
    end
end

