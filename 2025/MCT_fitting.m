clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 코드1에서 저장한 결과 불러오기
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

% Q_batt 후보 (예: 50 ~ 60 사이 10개 구간 대신, 여기서는 단일 값 사용)
Q_batt_candidates = 56.2396;  % [Ah]라고 가정

%% 3) fmincon 옵션 및 파라미터 초기값/제약조건 설정
%  5개 파라미터: X = [R0, R1, R2, tau1, tau2]
x0 = [0.001,  ... R0
      0.0005, ... R1
      0.0005, ... R2
      6.04,   ... tau1
      65];    ... tau2

lb = [0, 0, 0, 0, 10];
ub = [inf, inf, inf, 18.12, 195];

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
idx_firstNonZero = find(cellCurrent ~= 0, 1, 'first');

if isempty(idx_firstNonZero)
    % 전 구간 전류가 0이면 맨 처음(또는 맨 끝)을 택해도 되지만,
    % 여기서는 맨 앞 인덱스로 설정
    idx_init = 1;
    warning('cellCurrent가 전체 구간에서 0입니다. 초기 인덱스를 1로 설정했습니다.');
else
    % '처음으로 전류가 0이 아닌' 바로 이전 인덱스가 초기 구간의 마지막
    idx_init = idx_firstNonZero - 1;
    if idx_init < 1
        % 혹시 첫 샘플부터 전류가 0이 아닌 경우 예외처리
        idx_init = 1;
        warning('전류가 첫 지점부터 0이 아니므로, 초기 인덱스를 1로 설정했습니다.');
    end
end

% 초기 전압 및 초기 SOC
cellVoltage_init = cellVoltage_meas(idx_init);
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear', 'extrap');

% 시계열 간격
dt = [0; diff(time_s)];

%% 5) Q_batt 후보별로 (R0, R1, R2, tau1, tau2) 최적화
rmseArray = zeros(length(Q_batt_candidates), 1);
paramArray = zeros(length(Q_batt_candidates), 5);  % [R0, R1, R2, tau1, tau2]

for iQ = 1:length(Q_batt_candidates)
    Qbatt = Q_batt_candidates(iQ);  % [Ah]
    
    % (A) SOC(t) 계산
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

%% 6.1) FittingResult Table 생성
FittingResult = table(Q_batt_candidates(:), rmseArray, ...
                     paramArray(:, 1), paramArray(:, 2), paramArray(:, 3), ...
                     paramArray(:, 4), paramArray(:, 5), ...
   'VariableNames', {'Q_batt', 'RMSE', 'R0', 'R1', 'R2', 'tau1', 'tau2'});
disp('--- FittingResult Table ---');
disp(FittingResult);

%% 7) Q_batt vs RMSE 플롯
figure('Name', ['MCT-' num2str(mctNumber) ' : RMSE vs Q_batt (2-RC)'], ...
       'NumberTitle', 'off');
plot(Q_batt_candidates, rmseArray, 'o--', 'LineWidth', 1.2); hold on;
plot(bestQbatt, bestRMSE, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Q_{batt} [Ah]');
ylabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : 2-RC model RMSE vs Capacity']);
grid on;

%% 8) 최적 파라미터로 전압 & 전류 비교 플롯 (yyaxis)
% SOC 재계산
charge_integral_best = cumtrapz(time_s, cellCurrent);
SOC_t_best = SOC0 - (charge_integral_best / (bestQbatt * 3600)) * 100;

% 모델 전압 계산
V_estBest = modelVoltage_2RC_tau(bestParams, SOC_t_best, cellCurrent, dt, socOCV, ocvCellVoltage);

% 색상 정의
c = lines(3);

figure('Name', ['MCT-' num2str(mctNumber) ' : Voltage & Current (2-RC)'], ...
       'NumberTitle','off');

% (상단 Subplot) 전압(좌측) + 전류(우측)
subplot(2,1,1);
yyaxis left
plot(time_s, cellVoltage_meas, 'Color', c(2,:), 'LineWidth', 1.2); hold on;
plot(time_s, V_estBest,        'Color', c(3,:), 'LineWidth', 1.2);
ylabel('Cell Voltage (V)');

% 왼쪽 축 색상: 빨간색
ax = gca;
ax.YColor = 'r';

yyaxis right
plot(time_s, cellCurrent, 'Color', c(1,:), 'LineWidth', 1.2);
ylabel('Cell Current (A)');

% 오른쪽 축 색상: 파란색
ax = gca;
ax.YColor = 'b';

legend('V_{data}','V_{model}','I_{data}','Location','best');
xlabel('Time (s)');
% title(['MCT-' num2str(mctNumber) ...
%       ' : Voltage & Current (Q_{batt} = ' num2str(bestQbatt,'%.2f') 'Ah)']);
grid on;

% (하단 Subplot) 전압 오차
subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_estBest, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
% title('Error = (Measured - Model)');
grid on;

%% 9) 30초 간격으로 잘라서 여러 개 플롯 (Zoom In)
interval30 = 30;  % 30초 간격
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
    
    % 전압+전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
    plot(time_s(idx), V_estBest(idx),        'Color', c_local(3,:), 'LineWidth', 1.2);
    ylabel('Cell Voltage (V)');

    % 왼쪽 축 색상: 빨간색
    ax = gca;
    ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx),      'Color', c_local(1,:), 'LineWidth', 1.2);
    ylabel('Cell Current (A)');

    % 오른쪽 축 색상: 파란색
    ax = gca;
    ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    % title(['전압+전류 (30s 구간 #' num2str(iWin) ')']);
    grid on;
    
    % 에러
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_estBest(idx), 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    % title('Error = (Measured - Model)');
    grid on;
end

%% 10) 100초 간격으로 잘라서 여러 개 플롯 (Zoom In)
interval300 = 100;  % 100초 간격
numWindows300 = floor(maxTime / interval300);

for iWin = 1:numWindows300
    tStart = (iWin-1)*interval300;
    tEnd   = iWin*interval300;
    
    idx = (time_s >= tStart) & (time_s < tEnd);
    
    figure('Name', ['[100s 구간] #' num2str(iWin) ...
        ' (' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
        'NumberTitle', 'off');
    
    c_local = lines(3);

    % 전압+전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
    plot(time_s(idx), V_estBest(idx),        'Color', c_local(3,:), 'LineWidth', 1.2);
    ylabel('Cell Voltage (V)');

    % 왼쪽 축 색상: 빨간색
    ax = gca;
    ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx),      'Color', c_local(1,:), 'LineWidth', 1.2);
    ylabel('Cell Current (A)');

    % 오른쪽 축 색상: 파란색
    ax = gca;
    ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    % title(['전압+전류 (100s 구간 #' num2str(iWin) ')']);
    grid on;
    
    % 에러
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_estBest(idx), 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    % title('Error = (Measured - Model)');
    grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (아래는 로컬 함수 정의부) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    % 2개의 RC 전압 항 (초기값)
    Vrc1 = 0;
    Vrc2 = 0;
    
    N = length(SOC_t);
    V_est = zeros(N, 1);
    
    for k = 1:N
        % (1) OCV (SOC에 따른 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');
        
        % (2) R0 전압 강하
        IR0 = R0 * I_cell(k);
        
        % (3) RC1, RC2 업데이트
        if k > 1
            alpha1 = exp(-dt(k) / tau1);
            alpha2 = exp(-dt(k) / tau2);
            
            Vrc1 = Vrc1 * alpha1 + R1 * (1 - alpha1) * I_cell(k);
            Vrc2 = Vrc2 * alpha2 + R2 * (1 - alpha2) * I_cell(k);
        end
        
        % (4) 최종 모델 전압
        V_est(k) = OCV_now - IR0 - Vrc1 - Vrc2;
    end
end

