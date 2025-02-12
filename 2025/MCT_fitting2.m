%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [전체 코드 예시 - 1 RC 버전 (전류표시, 색상지정, TITLE 주석, 30s/100s 줌인)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 코드1에서 저장한 결과 불러오기
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

% Q_batt 후보 (여기서는 단일 값 가정)
Q_batt_candidates = 56.2396;  % 예: linspace(50,60,10);

%% 3) fmincon 옵션 및 파라미터 초기값/제약조건 설정
%  3개 파라미터: X = [R0, R1, tau1]
x0 = [0.001,  ... R0 초기값
      0.0005, ... R1 초기값
      180];    ... tau1 초기값

lb = [0, 0, 10];
ub = [inf, inf, 180];

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

%% 5) Q_batt 후보별로 (R0, R1, tau1) 최적화
rmseArray = zeros(length(Q_batt_candidates), 1);
paramArray = zeros(length(Q_batt_candidates), 3);  % [R0, R1, tau1]

for iQ = 1:length(Q_batt_candidates)
    Qbatt = Q_batt_candidates(iQ);  % [Ah]
    
    % (A) SOC(t) 계산 (+면 방전)
    charge_integral = cumtrapz(time_s, cellCurrent); % [A·s]
    SOC_t = SOC0 - (charge_integral / (Qbatt * 3600)) * 100;
    
    % (B) fmincon으로 최적화 (RMSE 최소화)
    costFunc = @(X) computeRMSE_1RC_tau(X, SOC_t, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
    [xOpt, fVal] = fmincon(costFunc, x0, [], [], [], [], lb, ub, [], options);
    
    rmseArray(iQ)    = fVal;
    paramArray(iQ,:) = xOpt;  % (R0, R1, tau1)
end

% (C) 최소 RMSE 지점
[bestRMSE, idxMin] = min(rmseArray);
bestQbatt = Q_batt_candidates(idxMin);
bestParams = paramArray(idxMin, :);  % [R0, R1, tau1]
bestR0   = bestParams(1);
bestR1   = bestParams(2);
bestTau1 = bestParams(3);

%% 6) 결과 출력
fprintf('\n=== MCT-%d 결과 (1-RC, 셀 전압; R0,R1,tau1) ===\n', mctNumber);
fprintf('Q_batt 후보: '); fprintf('%.1f ', Q_batt_candidates); fprintf('\n');
fprintf('RMSE들   : '); fprintf('%.4f ', rmseArray); fprintf('\n');
fprintf('--> 최소 RMSE = %.4f @ Q_batt = %.2f [Ah]\n', bestRMSE, bestQbatt);
fprintf('    (R0 = %.4f, R1 = %.4f, tau1 = %.1f)\n', bestR0, bestR1, bestTau1);

%% 6.1) FittingResult Table
FittingResult = table(Q_batt_candidates(:), rmseArray, ...
                     paramArray(:, 1), paramArray(:, 2), paramArray(:, 3), ...
                     'VariableNames', {'Q_batt','RMSE','R0','R1','tau1'});
disp('--- FittingResult Table ---');
disp(FittingResult);

%% 7) Q_batt vs RMSE 플롯
figure('Name', ['MCT-' num2str(mctNumber) ' : RMSE vs Q_batt (1-RC)'], ...
       'NumberTitle', 'off');
plot(Q_batt_candidates, rmseArray, 'o--', 'LineWidth', 1.2); hold on;
plot(bestQbatt, bestRMSE, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Q_{batt} [Ah]');
ylabel('RMSE (V)');
% title(['MCT-' num2str(mctNumber) ' : 1-RC 모델 RMSE vs Capacity']);  % <-- 주석처리
grid on;

%% 8) 최적 파라미터로 전압 & 전류 비교 플롯 (yyaxis)
% SOC 재계산
charge_integral_best = cumtrapz(time_s, cellCurrent);
SOC_t_best = SOC0 - (charge_integral_best / (bestQbatt * 3600)) * 100;

% 모델 전압 계산 (1-RC)
V_estBest = modelVoltage_1RC_tau(bestParams, SOC_t_best, cellCurrent, dt, ...
                                 socOCV, ocvCellVoltage);

% 색상 설정
c = lines(3);

figure('Name', ['MCT-' num2str(mctNumber) ' : Voltage & Current (1-RC)'], ...
       'NumberTitle','off');

% (상단 Subplot) 전압(좌) + 전류(우)
subplot(2,1,1);
yyaxis left
plot(time_s, cellVoltage_meas, 'Color', c(2,:), 'LineWidth', 1.2); hold on;
plot(time_s, V_estBest,        'Color', c(3,:), 'LineWidth', 1.2);
ylabel('Cell Voltage (V)');

yyaxis right
plot(time_s, cellCurrent,      'Color', c(1,:), 'LineWidth', 1.2);
ylabel('Cell Current (A)');

legend('V_{data}','V_{model}','I_{data}','Location','best');
xlabel('Time (s)');
% title(['MCT-' num2str(mctNumber) ...
%       ' 전압 비교 (1-RC, Q_{batt} = ' num2str(bestQbatt,'%.2f') 'Ah)']); % <-- 주석처리
grid on;

% (하단 Subplot) 전압 오차
subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_estBest, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
% title('Error = (Measured - Model)');  % <-- 주석처리
grid on;

%% 9) 30초 간격 줌인
interval30 = 30;  
maxTime = max(time_s);
numWindows30 = floor(maxTime / interval30);

for iWin = 1:numWindows30
    tStart = (iWin-1)*interval30;
    tEnd   = iWin*interval30;
    idx = (time_s >= tStart) & (time_s < tEnd);
    
    figure('Name', ['[30s 구간] #' num2str(iWin) ...
        ' (' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
        'NumberTitle', 'off');
    
    c_local = lines(3);

    % 전압 & 전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
    plot(time_s(idx), V_estBest(idx),        'Color', c_local(3,:), 'LineWidth', 1.2);
    ylabel('Cell Voltage (V)');

    yyaxis right
    plot(time_s(idx), cellCurrent(idx),      'Color', c_local(1,:), 'LineWidth', 1.2);
    ylabel('Cell Current (A)');

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    % title(['전압+전류 (30s 구간 #' num2str(iWin) ')']);  % <-- 주석처리
    grid on;

    % 에러
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_estBest(idx), 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    % title('Error = (Measured - Model)');  % <-- 주석처리
    grid on;
end

%% 10) 100초 간격 줌인
interval300 = 100;  % 요청에 따라 100초라면 바꾸면 됨
numWindows300 = floor(maxTime / interval300);

for iWin = 1:numWindows300
    tStart = (iWin-1)*interval300;
    tEnd   = iWin*interval300;
    idx = (time_s >= tStart) & (time_s < tEnd);
    
    figure('Name', ['[100s 구간] #' num2str(iWin) ...
        ' (' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
        'NumberTitle', 'off');
    
    c_local = lines(3);

    % 전압 & 전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
    plot(time_s(idx), V_estBest(idx),        'Color', c_local(3,:), 'LineWidth', 1.2);
    ylabel('Cell Voltage (V)');

    yyaxis right
    plot(time_s(idx), cellCurrent(idx),      'Color', c_local(1,:), 'LineWidth', 1.2);
    ylabel('Cell Current (A)');

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    % title(['전압+전류 (100s 구간 #' num2str(iWin) ')']);  % <-- 주석처리
    grid on;

    % 에러
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_estBest(idx), 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    % title('Error = (Measured - Model)');  % <-- 주석처리
    grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (아래는 로컬 함수 정의부) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = computeRMSE_1RC_tau(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    % X: [R0, R1, tau1]
    V_est = modelVoltage_1RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost = sqrt(mean((V_est - V_meas).^2));
end

function V_est = modelVoltage_1RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % X: [R0, R1, tau1]
    R0   = X(1);
    R1   = X(2);
    tau1 = X(3);
    
    % 1개의 RC 전압 항 (초기값)
    Vrc1 = 0;
    
    N = length(SOC_t);
    V_est = zeros(N, 1);
    
    for k = 1:N
        % (1) OCV (SOC에 따른 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');
        
        % (2) R0 전압 강하
        IR0 = R0 * I_cell(k);
        
        % (3) RC1 요소 업데이트
        if k > 1
            alpha1 = exp(-dt(k) / tau1);
            Vrc1 = Vrc1 * alpha1 + R1*(1 - alpha1)*I_cell(k);
        end
        
        % (4) 최종 모델 전압
        V_est(k) = OCV_now - IR0 - Vrc1;
    end
end

