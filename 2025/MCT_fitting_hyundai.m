%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [전체 코드] CC SOC vs BMS SOC (2-RC 모델 피팅)
%% + SOC 비교그래프 + 전체 시계열 2x2 플롯 + [30초/100초] 구간별 2x2 플롯
%% (왼쪽: CC SOC, 오른쪽: BMS SOC, 아래: Error 플롯)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) 사용자 입력: MCT 번호
mctNumber = input('피팅할 MCT 번호를 입력하세요 (예: 1~6): ');

%% 2) 데이터 로드 (예: MCT_Results.mat)
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

% (예시) 배터리 용량 [Ah]
Q_batt = 56.2396;  

%% 3) fmincon 옵션 및 초기값/제약조건 설정
% 2-RC: 파라미터 X = [R0, R1, R2, tau1, tau2]
x0 = [0.001,  ... R0
      0.0005, ... R1
      0.0005, ... R2
      6.04,   ... tau1
      65];    ... tau2

lb = [0,    0,    0,    0,    10 ];
ub = [inf, inf, inf, 18.12, 195];

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxIterations',1000,'MaxFunctionEvaluations',5000);

%% 4) MCT 데이터 불러오기
dataMCT = mctCellData{mctNumber};

time_s      = dataMCT.Time_s;
packVoltage = dataMCT.PackVoltage_V;
packCurrent = dataMCT.Current_A;   % (양수: 방전)
socDecimal  = dataMCT.SOC_decimal; 
socInteger  = dataMCT.SOC_integer; 

% (A) BMS SOC
SOC_bms = socInteger + socDecimal;  % [%]

% (B) Pack -> Cell 전압 변환 (192 직렬)
cellVoltage_meas = packVoltage / 192;
% (C) Pack -> Cell 전류 변환 (2 병렬)
cellCurrent = packCurrent / 2;

% (D) 시간 간격, 전체 시간
dt = [0; diff(time_s)];
maxTime = max(time_s);

%% 5) (A) 쿨롱카운팅 SOC 계산
% 5.1) 초기 SOC (OCV 기반)
idx_firstNonZero = find(cellCurrent ~= 0, 1, 'first');
if isempty(idx_firstNonZero)
    idx_init = 1;
    warning('전류가 전체 구간에서 0입니다. 초기 인덱스를 1로 설정합니다.');
else
    idx_init = idx_firstNonZero - 1;
    if idx_init < 1
        idx_init = 1;
        warning('첫 지점부터 전류가 0이 아니므로, 초기 인덱스를 1로 설정합니다.');
    end
end

cellVoltage_init = cellVoltage_meas(idx_init);
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear','extrap');

% 5.2) 누적 전류 적분 → SOC(t)
charge_integral = cumtrapz(time_s, cellCurrent);  % [A·s]
SOC_cc = SOC0 - (charge_integral/(Q_batt*3600))*100;  % [%]

%% 6) 2-RC 모델 파라미터 피팅
% (A) CC SOC 기반
costFunc_cc = @(X) computeRMSE_2RC_tau(X, SOC_cc, cellCurrent, dt, ...
                                       cellVoltage_meas, socOCV, ocvCellVoltage);
[bestParams_cc, bestRMSE_cc] = fmincon(costFunc_cc, x0, [], [], [], [], lb, ub, [], options);
bestQbatt_cc = Q_batt;

% (B) BMS SOC 기반
costFunc_bms = @(X) computeRMSE_2RC_tau(X, SOC_bms, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
[bestParams_bms, bestRMSE_bms] = fmincon(costFunc_bms, x0, [], [], [], [], lb, ub, [], options);
bestQbatt_bms = Q_batt;

%% 7) 결과 콘솔 출력
fprintf('\n=== MCT-%d 결과 비교 ===\n', mctNumber);

fprintf('--- CC SOC ---\n');
fprintf('RMSE = %.4f (V)\n', bestRMSE_cc);
fprintf('R0=%.5f, R1=%.5f, R2=%.5f, tau1=%.3f, tau2=%.3f\n', ...
    bestParams_cc(1), bestParams_cc(2), bestParams_cc(3), ...
    bestParams_cc(4), bestParams_cc(5));

fprintf('\n--- BMS SOC ---\n');
fprintf('RMSE = %.4f (V)\n', bestRMSE_bms);
fprintf('R0=%.5f, R1=%.5f, R2=%.5f, tau1=%.3f, tau2=%.3f\n', ...
    bestParams_bms(1), bestParams_bms(2), bestParams_bms(3), ...
    bestParams_bms(4), bestParams_bms(5));

%% 8) 결과 테이블 (2줄)
Method = {'CC SOC'; 'BMS SOC'};
FittingResult = table( ...
    Method, ...
    [bestQbatt_cc; bestQbatt_bms], ...
    [bestRMSE_cc;  bestRMSE_bms], ...
    [bestParams_cc(1); bestParams_bms(1)], ...
    [bestParams_cc(2); bestParams_bms(2)], ...
    [bestParams_cc(3); bestParams_bms(3)], ...
    [bestParams_cc(4); bestParams_bms(4)], ...
    [bestParams_cc(5); bestParams_bms(5)], ...
    'VariableNames', {'Method','Q_batt','RMSE','R0','R1','R2','tau1','tau2'} ...
);

disp('=== FittingResult Table (CC vs BMS) ===');
disp(FittingResult);

%% 9) SOC 비교 플롯 (CC SOC vs BMS SOC)
figure('Name','[SOC 비교] CC vs BMS','NumberTitle','off');
plot(time_s, SOC_cc,  'LineWidth',1.2); hold on;
plot(time_s, SOC_bms, 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('SOC (%)');
legend('CC SOC','BMS SOC','Location','best');
title(['MCT-' num2str(mctNumber) ' : CC SOC vs BMS SOC']);
grid on;

%% 10) 전체 시계열 피팅 결과 (2×2 Subplot, 왼=CC, 오른=BMS)
% 모델 전압 계산
V_est_cc  = modelVoltage_2RC_tau(bestParams_cc,  SOC_cc,  cellCurrent, dt, socOCV, ocvCellVoltage);
V_est_bms = modelVoltage_2RC_tau(bestParams_bms, SOC_bms, cellCurrent, dt, socOCV, ocvCellVoltage);

figure('Name','[Fitting - 전체 2x2] CC vs BMS','NumberTitle','off');

% (1,1) CC SOC: V_meas vs V_model
subplot(2,2,1);
plot(time_s, cellVoltage_meas, 'b-', 'LineWidth',1.1); hold on;
plot(time_s, V_est_cc,        'r--','LineWidth',1.1);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('CC SOC: V_{meas} vs V_{model}');
legend('V_{meas}','V_{model}','Location','best');
grid on;

% (1,2) BMS SOC: V_meas vs V_model
subplot(2,2,2);
plot(time_s, cellVoltage_meas, 'b-', 'LineWidth',1.1); hold on;
plot(time_s, V_est_bms,       'r--','LineWidth',1.1);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('BMS SOC: V_{meas} vs V_{model}');
legend('V_{meas}','V_{model}','Location','best');
grid on;

% (2,1) CC SOC: Error
subplot(2,2,3);
plot(time_s, cellVoltage_meas - V_est_cc, 'k-','LineWidth',1.1);
xlabel('Time (s)'); ylabel('Voltage Error (V)');
title('CC SOC: Error');
grid on;

% (2,2) BMS SOC: Error
subplot(2,2,4);
plot(time_s, cellVoltage_meas - V_est_bms, 'k-','LineWidth',1.1);
xlabel('Time (s)'); ylabel('Voltage Error (V)');
title('BMS SOC: Error');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11) [30초 간격] 구간별 2×2 Subplot (왼=CC, 오른=BMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interval30 = 30;
numWindows30 = floor(maxTime / interval30);

for iWin = 1:numWindows30
    tStart = (iWin-1)*interval30;
    tEnd   = iWin*interval30;
    idx = (time_s >= tStart) & (time_s < tEnd);
    
    figure('Name', ...
        ['[30s 구간 #' num2str(iWin) '] MCT-' num2str(mctNumber)], ...
        'NumberTitle','off');
    
    % (1,1) CC SOC: V_meas vs V_model
    subplot(2,2,1);
    plot(time_s(idx), cellVoltage_meas(idx), 'b-','LineWidth',1.1); hold on;
    plot(time_s(idx), V_est_cc(idx),        'r--','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Voltage (V)');
    title(['CC SOC: ' num2str(tStart) '~' num2str(tEnd) ' s']);
    legend('V_{meas}','V_{model}','Location','best');
    grid on;

    % (1,2) BMS SOC: V_meas vs V_model
    subplot(2,2,2);
    plot(time_s(idx), cellVoltage_meas(idx), 'b-','LineWidth',1.1); hold on;
    plot(time_s(idx), V_est_bms(idx),       'r--','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Voltage (V)');
    title(['BMS SOC: ' num2str(tStart) '~' num2str(tEnd) ' s']);
    legend('V_{meas}','V_{model}','Location','best');
    grid on;

    % (2,1) CC SOC: Error
    subplot(2,2,3);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_cc(idx), 'k-','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Error (V)');
    title('CC SOC: Error');
    grid on;

    % (2,2) BMS SOC: Error
    subplot(2,2,4);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_bms(idx), 'k-','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Error (V)');
    title('BMS SOC: Error');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 12) [100초 간격] 구간별 2×2 Subplot (왼=CC, 오른=BMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interval100 = 100;
numWindows100 = floor(maxTime / interval100);

for iWin = 1:numWindows100
    tStart = (iWin-1)*interval100;
    tEnd   = iWin*interval100;
    idx = (time_s >= tStart) & (time_s < tEnd);
    
    figure('Name', ...
        ['[100s 구간 #' num2str(iWin) '] MCT-' num2str(mctNumber)], ...
        'NumberTitle','off');
    
    % (1,1) CC SOC: V_meas vs V_model
    subplot(2,2,1);
    plot(time_s(idx), cellVoltage_meas(idx), 'b-','LineWidth',1.1); hold on;
    plot(time_s(idx), V_est_cc(idx),        'r--','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Voltage (V)');
    title(['CC SOC: ' num2str(tStart) '~' num2str(tEnd) ' s']);
    legend('V_{meas}','V_{model}','Location','best');
    grid on;

    % (1,2) BMS SOC: V_meas vs V_model
    subplot(2,2,2);
    plot(time_s(idx), cellVoltage_meas(idx), 'b-','LineWidth',1.1); hold on;
    plot(time_s(idx), V_est_bms(idx),       'r--','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Voltage (V)');
    title(['BMS SOC: ' num2str(tStart) '~' num2str(tEnd) ' s']);
    legend('V_{meas}','V_{model}','Location','best');
    grid on;

    % (2,1) CC SOC: Error
    subplot(2,2,3);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_cc(idx), 'k-','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Error (V)');
    title('CC SOC: Error');
    grid on;

    % (2,2) BMS SOC: Error
    subplot(2,2,4);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_bms(idx), 'k-','LineWidth',1.1);
    xlabel('Time (s)'); ylabel('Error (V)');
    title('BMS SOC: Error');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 로컬 함수 정의
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = computeRMSE_2RC_tau(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    % X: [R0, R1, R2, tau1, tau2]
    % RMSE = sqrt( mean( (V_meas - V_est)^2 ) )
    
    V_est = modelVoltage_2RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost  = sqrt(mean((V_meas - V_est).^2));
end

function V_est = modelVoltage_2RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % 2-RC 모델 식:
    %  V_est(k) = OCV(SOC(k)) - R0*I(k) - Vrc1(k) - Vrc2(k)
    %  Vrc1(k+1) = exp(-dt/tau1)*Vrc1(k) + R1*(1-exp(-dt/tau1))*I(k)
    %  Vrc2(k+1) = exp(-dt/tau2)*Vrc2(k) + R2*(1-exp(-dt/tau2))*I(k)
    
    R0   = X(1);
    R1   = X(2);
    R2   = X(3);
    tau1 = X(4);
    tau2 = X(5);
    
    N = length(SOC_t);
    V_est = zeros(N,1);
    
    Vrc1 = 0;  % RC1 초기값
    Vrc2 = 0;  % RC2 초기값
    
    for k = 1:N
        % (1) OCV (SOC->전압 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear','extrap');
        
        % (2) R0 전압강하
        IR0 = R0 * I_cell(k);
        
        % (3) RC 항 업데이트
        if k > 1
            alpha1 = exp(-dt(k)/tau1);
            alpha2 = exp(-dt(k)/tau2);
            Vrc1 = alpha1 * Vrc1 + R1*(1 - alpha1)*I_cell(k);
            Vrc2 = alpha2 * Vrc2 + R2*(1 - alpha2)*I_cell(k);
        end
        
        % (4) 최종 모델 전압
        V_est(k) = OCV_now - IR0 - Vrc1 - Vrc2;
    end
end



