%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 코드2 : 특정 MCT 번호만 입력받아,
%         2-RC 모델 (R0 + R1||C1 + R2||C2) 파라미터 (5개) 피팅 + 1-RC 모델 추가
%         - packVoltage -> cellVoltage (/192)
%         - packCurrent -> cellCurrent (/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

if ~(mctNumber >= 1 && mctNumber <= 6)
    error('MCT 번호는 1부터 6 사이의 정수여야 합니다.');
end

%% 2) 코드1에서 저장한 결과 불러오기 (예시)
%    실제로는 "MCT_Results.mat" 파일에 맞게 경로/변수명을 수정해야 합니다.
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

% Q_batt 후보 (예: 50 ~ 60 사이 10개 구간)
Q_batt_candidates = linspace(50,60,10);  % [Ah]라고 가정

%% 3) fmincon 옵션 및 파라미터 초기값/제약조건 설정

%----------- (A) 2-RC 모델용 -----------
%  5개 파라미터: X_2RC = [R0, R1, C1, R2, C2]
lb_2RC = [0, 0,    0,   0,    0];
ub_2RC = [inf, inf, inf, inf, inf];
x0_2RC = [0.038, 0.02, 150, 0.02, 1000];  % 초기 추정값 

%----------- (B) 1-RC 모델용 -----------
%  3개 파라미터: X_1RC = [R0, R1, C1]
lb_1RC = [0, 0, 0];
ub_1RC = [inf, inf, inf];
x0_1RC = [0.014, 0.08, 1500];  % 초기 추정값(적절히 조정)

%----------- fmincon 옵션 -----------
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
%     예: 첫 5개 샘플 평균 전압으로 SOC0 추정
cellVoltage_init = mean(cellVoltage_meas(1:5));  % 실제 데이터에 맞게 조정 가능

% 만약 ocvCellVoltage가 단조 정렬되어 있다면 그대로 사용
% 단조 정렬이 불분명하면 uCellVoltage, uSocOCV(정렬+unique) 사용
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear', 'extrap');

% 시계열 간격
dt = [0; diff(time_s)];

%% =========================================================================
%              (I) 2-RC 모델 피팅 (기존 코드와 동일)
%% =========================================================================

rmseArray_2RC = zeros(length(Q_batt_candidates), 1);
paramArray_2RC = zeros(length(Q_batt_candidates), 5);  % [R0, R1, C1, R2, C2]

for iQ = 1:length(Q_batt_candidates)
    Qbatt = Q_batt_candidates(iQ);  % [Ah]
    
    % (A) SOC(t) 계산
    charge_integral = cumtrapz(time_s, cellCurrent); % [A·s]
    SOC_t = SOC0 - (charge_integral / (Qbatt * 3600)) * 100;
    
    % (B) fmincon으로 최적화 (RMSE 최소화)
    costFunc_2RC = @(X) computeRMSE_2RC(X, SOC_t, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
    [xOpt, fVal] = fmincon(costFunc_2RC, x0_2RC, [], [], [], [], lb_2RC, ub_2RC, [], options);
    
    rmseArray_2RC(iQ) = fVal;
    paramArray_2RC(iQ, :) = xOpt;  % (R0, R1, C1, R2, C2)
end

% (C) 최소 RMSE 지점
[bestRMSE_2RC, idxMin_2RC] = min(rmseArray_2RC);
bestQbatt_2RC = Q_batt_candidates(idxMin_2RC);
bestParams_2RC = paramArray_2RC(idxMin_2RC, :);  % 5개 파라미터
bestR0_2RC = bestParams_2RC(1);
bestR1_2RC = bestParams_2RC(2);
bestC1_2RC = bestParams_2RC(3);
bestR2_2RC = bestParams_2RC(4);
bestC2_2RC = bestParams_2RC(5);

%% 결과 출력 (2-RC)
fprintf('\n=== MCT-%d 결과 (2-RC, 셀 전압) ===\n', mctNumber);
fprintf('Q_batt 후보: '); fprintf('%.1f ', Q_batt_candidates); fprintf('\n');
fprintf('RMSE들   : '); fprintf('%.4f ', rmseArray_2RC); fprintf('\n');
fprintf('--> 최소 RMSE = %.4f @ Q_batt = %.2f [Ah]\n', bestRMSE_2RC, bestQbatt_2RC);
fprintf('    (R0 = %.4f, R1 = %.4f, C1 = %.1f, R2 = %.4f, C2 = %.1f)\n', ...
        bestR0_2RC, bestR1_2RC, bestC1_2RC, bestR2_2RC, bestC2_2RC);

%% 결과 Table (2-RC)
FittingResult_2RC = table(Q_batt_candidates(:), rmseArray_2RC, ...
                          paramArray_2RC(:, 1), paramArray_2RC(:, 2), paramArray_2RC(:, 3), ...
                          paramArray_2RC(:, 4), paramArray_2RC(:, 5), ...
                         'VariableNames', {'Q_batt', 'RMSE', 'R0', 'R1', 'C1', 'R2', 'C2'});

disp('--- FittingResult Table (2-RC) ---');
disp(FittingResult_2RC);

%% RMSE vs Q_batt 플롯 (2-RC)
figure('Name', ['MCT-' num2str(mctNumber) ' : RMSE vs Q_batt (2-RC)'], 'NumberTitle', 'off');
plot(Q_batt_candidates, rmseArray_2RC, 'o--', 'LineWidth', 1.2);
hold on;
plot(bestQbatt_2RC, bestRMSE_2RC, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Q_{batt} [Ah]');
ylabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : 2-RC 모델 RMSE vs Capacity']);
grid on;

%% 최적 파라미터로 전압 비교 플롯 (2-RC)
charge_integral_best_2RC = cumtrapz(time_s, cellCurrent);
SOC_t_best_2RC = SOC0 - (charge_integral_best_2RC / (bestQbatt_2RC * 3600)) * 100;
V_estBest_2RC = modelVoltage_2RC(bestParams_2RC, SOC_t_best_2RC, cellCurrent, dt, socOCV, ocvCellVoltage);

figure('Name', ['MCT-' num2str(mctNumber) ' : Voltage Comparison (2-RC)'], 'NumberTitle', 'off');
subplot(2,1,1);
plot(time_s, cellVoltage_meas, '-b', 'LineWidth', 1.2); hold on;
plot(time_s, V_estBest_2RC, '--r', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Cell Voltage (V)');
legend('Measured', 'Model', 'Location', 'best');
title(['MCT-' num2str(mctNumber) ' 전압 비교 (2-RC, Q_{batt} = ' num2str(bestQbatt_2RC, '%.2f') 'Ah)']);
grid on;

subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_estBest_2RC, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
title('Error = (Measured - Model)');
grid on;

%% =========================================================================
%              (II) 1-RC 모델 피팅 추가
%% =========================================================================

rmseArray_1RC = zeros(length(Q_batt_candidates), 1);
paramArray_1RC = zeros(length(Q_batt_candidates), 3);  % [R0, R1, C1]

for iQ = 1:length(Q_batt_candidates)
    Qbatt = Q_batt_candidates(iQ);  % [Ah]
    
    % (A) SOC(t) 계산
    charge_integral = cumtrapz(time_s, cellCurrent); % [A·s]
    SOC_t = SOC0 - (charge_integral / (Qbatt * 3600)) * 100;
    
    % (B) fmincon으로 최적화 (RMSE 최소화)
    costFunc_1RC = @(X) computeRMSE_1RC(X, SOC_t, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
    [xOpt, fVal] = fmincon(costFunc_1RC, x0_1RC, [], [], [], [], lb_1RC, ub_1RC, [], options);
    
    rmseArray_1RC(iQ) = fVal;
    paramArray_1RC(iQ, :) = xOpt;  % (R0, R1, C1)
end

% (C) 최소 RMSE 지점
[bestRMSE_1RC, idxMin_1RC] = min(rmseArray_1RC);
bestQbatt_1RC = Q_batt_candidates(idxMin_1RC);
bestParams_1RC = paramArray_1RC(idxMin_1RC, :);  % 3개 파라미터
bestR0_1RC = bestParams_1RC(1);
bestR1_1RC = bestParams_1RC(2);
bestC1_1RC = bestParams_1RC(3);

%% 결과 출력 (1-RC)
fprintf('\n=== MCT-%d 결과 (1-RC, 셀 전압) ===\n', mctNumber);
fprintf('Q_batt 후보: '); fprintf('%.1f ', Q_batt_candidates); fprintf('\n');
fprintf('RMSE들   : '); fprintf('%.4f ', rmseArray_1RC); fprintf('\n');
fprintf('--> 최소 RMSE = %.4f @ Q_batt = %.2f [Ah]\n', bestRMSE_1RC, bestQbatt_1RC);
fprintf('    (R0 = %.4f, R1 = %.4f, C1 = %.1f)\n', ...
        bestR0_1RC, bestR1_1RC, bestC1_1RC);

%% 결과 Table (1-RC)
FittingResult_1RC = table(Q_batt_candidates(:), rmseArray_1RC, ...
                          paramArray_1RC(:, 1), paramArray_1RC(:, 2), paramArray_1RC(:, 3), ...
                         'VariableNames', {'Q_batt', 'RMSE', 'R0', 'R1', 'C1'});
disp('--- FittingResult Table (1-RC) ---');
disp(FittingResult_1RC);

%% RMSE vs Q_batt 플롯 (1-RC)
figure('Name', ['MCT-' num2str(mctNumber) ' : RMSE vs Q_batt (1-RC)'], 'NumberTitle', 'off');
plot(Q_batt_candidates, rmseArray_1RC, 'o--', 'LineWidth', 1.2);
hold on;
plot(bestQbatt_1RC, bestRMSE_1RC, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Q_{batt} [Ah]');
ylabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : 1-RC 모델 RMSE vs Capacity']);
grid on;

%% 최적 파라미터로 전압 비교 플롯 (1-RC)
charge_integral_best_1RC = cumtrapz(time_s, cellCurrent);
SOC_t_best_1RC = SOC0 - (charge_integral_best_1RC / (bestQbatt_1RC * 3600)) * 100;
V_estBest_1RC = modelVoltage_1RC(bestParams_1RC, SOC_t_best_1RC, cellCurrent, dt, socOCV, ocvCellVoltage);

figure('Name', ['MCT-' num2str(mctNumber) ' : Voltage Comparison (1-RC)'], 'NumberTitle', 'off');
subplot(2,1,1);
plot(time_s, cellVoltage_meas, '-b', 'LineWidth', 1.2); hold on;
plot(time_s, V_estBest_1RC, '--r', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Cell Voltage (V)');
legend('Measured', 'Model', 'Location', 'best');
title(['MCT-' num2str(mctNumber) ' 전압 비교 (1-RC, Q_{batt} = ' num2str(bestQbatt_1RC, '%.2f') 'Ah)']);
grid on;

subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_estBest_1RC, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
title('Error = (Measured - Model)');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 서브펑션 (RMSE 계산 및 모델 전압 계산)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================== (1) 2-RC 서브펑션 ================================
function cost = computeRMSE_2RC(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    % X: [R0, R1, C1, R2, C2]
    V_est = modelVoltage_2RC(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost = sqrt(mean((V_est - V_meas).^2));
end

function V_est = modelVoltage_2RC(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % X = [R0, R1, C1, R2, C2]
    R0 = X(1);    R1 = X(2);    C1 = X(3);
    R2 = X(4);    C2 = X(5);
    
    % 2개의 RC 전압 항 (초기값)
    Vrc1 = 0;
    Vrc2 = 0;
    
    N = length(SOC_t);
    V_est = zeros(N, 1);
    
    for k = 1:N
        % (1) OCV 추출 (SOC에 따른 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');
        
        % (2) R0에 의한 전압 강하
        IR0 = R0 * I_cell(k);
        
        % (3) RC1, RC2 항 업데이트
        if k > 1
            alpha1 = exp(-dt(k) / (R1 * C1));
            alpha2 = exp(-dt(k) / (R2 * C2));
            Vrc1 = Vrc1 * alpha1 + R1 * (1 - alpha1) * I_cell(k);
            Vrc2 = Vrc2 * alpha2 + R2 * (1 - alpha2) * I_cell(k);
        end
        
        % (4) 최종 전압
        V_est(k) = OCV_now - IR0 - Vrc1 - Vrc2;
    end
end

%======================== (2) 1-RC 서브펑션 ================================
function cost = computeRMSE_1RC(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    % X: [R0, R1, C1]
    V_est = modelVoltage_1RC(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost = sqrt(mean((V_est - V_meas).^2));
end

function V_est = modelVoltage_1RC(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % X = [R0, R1, C1]
    R0 = X(1);    
    R1 = X(2);    
    C1 = X(3);
    
    % 1개의 RC 전압 항 (초기값)
    Vrc = 0;
    
    N = length(SOC_t);
    V_est = zeros(N, 1);
    
    for k = 1:N
        % (1) OCV 추출 (SOC에 따른 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');
        
        % (2) R0에 의한 전압 강하
        IR0 = R0 * I_cell(k);
        
        % (3) RC 항 업데이트
        if k > 1
            alpha = exp(-dt(k) / (R1 * C1));
            Vrc = Vrc * alpha + R1 * (1 - alpha) * I_cell(k);
        end
        
        % (4) 최종 전압
        V_est(k) = OCV_now - IR0 - Vrc;
    end
end
