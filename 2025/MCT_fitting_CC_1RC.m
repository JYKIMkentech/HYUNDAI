%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1RC 모델 전압 Fitting (R0, R1, tau1 추정) - 초기 RC 전압을 따로 선언 없이
% 첫 샘플에서 Vrc1 = R1*I(1)*(1 - exp(-dt(1)/tau1)) 로 바로 설정
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 데이터 불러오기
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

% Q_batt 후보 (단일 값 가정 예시)
Q_batt_candidates = 56.2396;  % 필요시 linspace(50,60,10) 등 사용 가능

%% 3) fmincon 옵션 및 파라미터 초기값/제약조건
x0 = [0.001,  ... R0 초기값
      0.0005, ... R1 초기값
      36];    ... tau1 초기값
lb = [0, 0, 0];
ub = [inf, inf, 108];

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                       'MaxIterations',1000,'MaxFunctionEvaluations',5000);

%% 4) MCT 데이터 준비
dataMCT = mctCellData{mctNumber};

time_s      = dataMCT.Time_s;
packVoltage = dataMCT.PackVoltage_V;
packCurrent = dataMCT.Current_A;   % + : 방전 전류

% Pack → Cell로 변환
cellVoltage_meas = packVoltage / 192;  % 직렬 192셀 가정
cellCurrent = packCurrent / 2;         % 병렬 2셀 가정

%% (중요) 초기 SOC 계산
idx_firstNonZero = find(cellCurrent ~= 0, 1, 'first');
if isempty(idx_firstNonZero)
    idx_init = 1;
    warning('cellCurrent가 전체 구간에서 0이므로 초기 인덱스=1로 설정');
else
    idx_init = idx_firstNonZero - 1;
    if idx_init < 1
        idx_init = 1;
        warning('전류가 첫 지점부터 0이 아니므로, 초기 인덱스=1로 설정');
    end
end

cellVoltage_init = cellVoltage_meas(idx_init);
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear', 'extrap');

dt = [1; diff(time_s)];

%% 5) Q_batt 후보별로 (R0, R1, tau1) 최적화
rmseArray = zeros(length(Q_batt_candidates), 1);
paramArray = zeros(length(Q_batt_candidates), 3);  % [R0, R1, tau1]

for iQ = 1:length(Q_batt_candidates)
    Qbatt = Q_batt_candidates(iQ);  % [Ah]
    
    % SOC(t) 계산 (+ 방전)
    charge_integral = cumtrapz(time_s, cellCurrent); % [A·s]
    SOC_t = SOC0 - (charge_integral / (Qbatt * 3600)) * 100;
    
    % fmincon으로 (R0,R1,tau1) 최적화
    costFunc = @(X) computeRMSE_1RC_tau(X, SOC_t, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
    [xOpt, fVal] = fmincon(costFunc, x0, [], [], [], [], lb, ub, [], options);
    
    rmseArray(iQ)    = fVal;
    paramArray(iQ,:) = xOpt;
end

% 최소 RMSE 지점
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

% 결과 Table
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
grid on;

%% 8) 최적 파라미터 전압 & 전류 비교
charge_integral_best = cumtrapz(time_s, cellCurrent);
SOC_t_best = SOC0 - (charge_integral_best / (bestQbatt * 3600)) * 100;

V_estBest = modelVoltage_1RC_tau(bestParams, SOC_t_best, cellCurrent, dt, ...
                                 socOCV, ocvCellVoltage);

c = lines(3);
figure('Name', ['MCT-' num2str(mctNumber) ' : Voltage & Current (1-RC)'], ...
       'NumberTitle','off');

subplot(2,1,1);
yyaxis left
plot(time_s, cellVoltage_meas, 'Color', c(2,:), 'LineWidth', 1.2); hold on;
plot(time_s, V_estBest,        'Color', c(3,:), 'LineWidth', 1.2);
ylabel('Cell Voltage (V)');
ax = gca;
ax.YColor = 'r';

yyaxis right
plot(time_s, cellCurrent,      'Color', c(1,:), 'LineWidth', 1.2);
ylabel('Cell Current (A)');
ax = gca;
ax.YColor = 'b';

legend('V_{data}','V_{model}','I_{data}','Location','best');
xlabel('Time (s)');
grid on;

subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_estBest, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
grid on;

%% 9) 30초 간격 줌인
interval30 = 30;  
maxTime = max(time_s);
numWindows30 = floor(maxTime / interval30);

for iWin = 1:numWindows30
    tStart = (iWin-1)*interval30;
    tEnd   = iWin*interval30;
    idx = (time_s >= tStart & time_s < tEnd);
    
    figure('Name', ['[30s 구간] #' num2str(iWin) ...
        ' (' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
        'NumberTitle', 'off');
    
    c_local = lines(3);
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
    plot(time_s(idx), V_estBest(idx),        'Color', c_local(3,:), 'LineWidth', 1.2);
    ylabel('Cell Voltage (V)');
    ax = gca;
    ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx),      'Color', c_local(1,:), 'LineWidth', 1.2);
    ylabel('Cell Current (A)');
    ax = gca;
    ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    grid on;

    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_estBest(idx), 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    grid on;
end

%% 10) 100초 간격 줌인
interval300 = 100;
numWindows300 = floor(maxTime / interval300);

for iWin = 1:numWindows300
    tStart = (iWin-1)*interval300;
    tEnd   = iWin*interval300;
    idx = (time_s >= tStart & time_s < tEnd);
    
    figure('Name', ['[100s 구간] #' num2str(iWin) ...
        ' (' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
        'NumberTitle', 'off');
    
    c_local = lines(3);
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
    plot(time_s(idx), V_estBest(idx),        'Color', c_local(3,:), 'LineWidth', 1.2);
    ylabel('Cell Voltage (V)');
    ax = gca;
    ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx),      'Color', c_local(1,:), 'LineWidth', 1.2);
    ylabel('Cell Current (A)');
    ax = gca;
    ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    grid on;

    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_estBest(idx), 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (로컬 함수들)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = computeRMSE_1RC_tau(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    V_est = modelVoltage_1RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost = sqrt(mean((V_est - V_meas).^2));
end

function V_est = modelVoltage_1RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % X: [R0, R1, tau1]
    R0   = X(1);
    R1   = X(2);
    tau1 = X(3);

    N = length(SOC_t);
    V_est = zeros(N, 1);

    for k = 1:N
        % (1) OCV
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');

        % (2) R0 전압강하
        IR0 = R0 * I_cell(k);

        % (3) RC1 업데이트
        alpha1 = exp(-dt(k)/tau1);
        if k == 1
            % 첫 샘플에서만 RC 전압을 직접 설정
            Vrc1 = R1 * I_cell(1) * (1 - alpha1);
        else
            % 이후는 기존 공식
            Vrc1 = Vrc1*alpha1 + R1*(1 - alpha1)*I_cell(k);
        end

        % (4) 최종 전압
        V_est(k) = OCV_now - IR0 - Vrc1;
    end
end
