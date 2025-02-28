%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-RC 모델 : BMS SOC vs CC SOC 비교
%%   - Subplot(2,1), lines(3) 색상 사용
%%   - 전체/30초/100초 간격
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
% 1-RC: 파라미터 X = [R0, R1, tau1]
x0 = [0.001,  ... R0
      0.0005, ... R1
      306];   ... tau1 (예시)

lb = [0, 0, 0];
ub = [inf, inf, 918];

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
dt = [1; diff(time_s)];  % 여기서는 첫 샘플 간격을 1초로 가정
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

%% 6) 1-RC 모델 파라미터 피팅 (BMS SOC / CC SOC)
% (A) CC SOC 기반
costFunc_cc = @(X) computeRMSE_1RC_tau(X, SOC_cc, cellCurrent, dt, ...
                                       cellVoltage_meas, socOCV, ocvCellVoltage);
[bestParams_cc, bestRMSE_cc] = fmincon(costFunc_cc, x0, [], [], [], [], lb, ub, [], options);

% (B) BMS SOC 기반
costFunc_bms = @(X) computeRMSE_1RC_tau(X, SOC_bms, cellCurrent, dt, ...
                                        cellVoltage_meas, socOCV, ocvCellVoltage);
[bestParams_bms, bestRMSE_bms] = fmincon(costFunc_bms, x0, [], [], [], [], lb, ub, [], options);

%% 7) 모델 전압 계산
V_est_cc  = modelVoltage_1RC_tau(bestParams_cc,  SOC_cc,  cellCurrent, dt, socOCV, ocvCellVoltage);
V_est_bms = modelVoltage_1RC_tau(bestParams_bms, SOC_bms, cellCurrent, dt, socOCV, ocvCellVoltage);

%% 8) 결과 콘솔 출력
fprintf('\n=== MCT-%d 결과 (1-RC 모델) ===\n', mctNumber);
fprintf('--- CC SOC ---\n');
fprintf('RMSE = %.4f (V)\n', bestRMSE_cc);
fprintf('R0=%.5f, R1=%.5f, tau1=%.3f\n', ...
    bestParams_cc(1), bestParams_cc(2), bestParams_cc(3));

fprintf('\n--- BMS SOC ---\n');
fprintf('RMSE = %.4f (V)\n', bestRMSE_bms);
fprintf('R0=%.5f, R1=%.5f, tau1=%.3f\n', ...
    bestParams_bms(1), bestParams_bms(2), bestParams_bms(3));

%% --- (선택) 결과 테이블 ---
Method = {'CC SOC'; 'BMS SOC'};
RMSE   = [bestRMSE_cc; bestRMSE_bms];
R0     = [bestParams_cc(1); bestParams_bms(1)];
R1     = [bestParams_cc(2); bestParams_bms(2)];
tau1   = [bestParams_cc(3); bestParams_bms(3)];
FittingResult = table(Method, RMSE, R0, R1, tau1);
disp(FittingResult);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [A] BMS SOC 기반 모델 : 전체 시간 플롯 (2×1 Subplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','[BMS SOC] Entire time','NumberTitle','off');
c_local = lines(3);

% (A-1) 전압 & 전류
subplot(2,1,1);
yyaxis left
plot(time_s, cellVoltage_meas, 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
plot(time_s, V_est_bms,       'Color', c_local(3,:), 'LineWidth', 1.2);
ylabel('Cell Voltage (V)');
ax = gca; ax.YColor = 'r';  % 왼쪽 축 빨간색으로

yyaxis right
plot(time_s, cellCurrent,     'Color', c_local(1,:), 'LineWidth', 1.2);
ylabel('Cell Current (A)');
ax = gca; ax.YColor = 'b';   % 오른쪽 축 파란색으로

legend('V_{data}','V_{model}','I_{data}','Location','best');
xlabel('Time (s)');
grid on;

% (A-2) 전압 Error
subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_est_bms, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [B] CC SOC 기반 모델 : 전체 시간 플롯 (2×1 Subplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','[CC SOC] Entire time','NumberTitle','off');
c_local = lines(3);

% (B-1) 전압 & 전류
subplot(2,1,1);
yyaxis left
plot(time_s, cellVoltage_meas, 'Color', c_local(2,:), 'LineWidth', 1.2); hold on;
plot(time_s, V_est_cc,        'Color', c_local(3,:), 'LineWidth', 1.2);
ylabel('Cell Voltage (V)');
ax = gca; ax.YColor = 'r';

yyaxis right
plot(time_s, cellCurrent,     'Color', c_local(1,:), 'LineWidth', 1.2);
ylabel('Cell Current (A)');
ax = gca; ax.YColor = 'b';

legend('V_{data}','V_{model}','I_{data}','Location','best');
xlabel('Time (s)');
grid on;

% (B-2) 전압 Error
subplot(2,1,2);
plot(time_s, cellVoltage_meas - V_est_cc, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage Error (V)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [C] BMS SOC 기반 모델 : 30초 간격 줌인 (2×1 Subplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval30 = 30;
numWindows30 = floor(maxTime / interval30);

for iWin = 1:numWindows30
    tStart = (iWin-1)*interval30;
    tEnd   = iWin*interval30;
    idx = (time_s >= tStart) & (time_s < tEnd);

    figure('Name',...
        ['[BMS, 30s 구간 #' num2str(iWin) '] ' ...
         '(' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
         'NumberTitle','off');
    c_local = lines(3);

    % (C-1) 전압 & 전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth',1.2); hold on;
    plot(time_s(idx), V_est_bms(idx),        'Color', c_local(3,:), 'LineWidth',1.2);
    ylabel('Cell Voltage (V)');
    ax = gca; ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx), 'Color', c_local(1,:), 'LineWidth',1.2);
    ylabel('Cell Current (A)');
    ax = gca; ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    grid on;

    % (C-2) 전압 Error
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_bms(idx), 'k', 'LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [D] BMS SOC 기반 모델 : 100초 간격 줌인 (2×1 Subplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval100 = 100;
numWindows100 = floor(maxTime / interval100);

for iWin = 1:numWindows100
    tStart = (iWin-1)*interval100;
    tEnd   = iWin*interval100;
    idx = (time_s >= tStart) & (time_s < tEnd);

    figure('Name',...
        ['[BMS, 100s 구간 #' num2str(iWin) '] ' ...
         '(' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
         'NumberTitle','off');
    c_local = lines(3);

    % (D-1) 전압 & 전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth',1.2); hold on;
    plot(time_s(idx), V_est_bms(idx),        'Color', c_local(3,:), 'LineWidth',1.2);
    ylabel('Cell Voltage (V)');
    ax = gca; ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx), 'Color', c_local(1,:), 'LineWidth',1.2);
    ylabel('Cell Current (A)');
    ax = gca; ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    grid on;

    % (D-2) 전압 Error
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_bms(idx), 'k', 'LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [E] CC SOC 기반 모델 : 30초 간격 줌인 (2×1 Subplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval30 = 30;
numWindows30 = floor(maxTime / interval30);

for iWin = 1:numWindows30
    tStart = (iWin-1)*interval30;
    tEnd   = iWin*interval30;
    idx = (time_s >= tStart) & (time_s < tEnd);

    figure('Name',...
        ['[CC, 30s 구간 #' num2str(iWin) '] ' ...
         '(' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
         'NumberTitle','off');
    c_local = lines(3);

    % (E-1) 전압 & 전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth',1.2); hold on;
    plot(time_s(idx), V_est_cc(idx),        'Color', c_local(3,:), 'LineWidth',1.2);
    ylabel('Cell Voltage (V)');
    ax = gca; ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx), 'Color', c_local(1,:), 'LineWidth',1.2);
    ylabel('Cell Current (A)');
    ax = gca; ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    grid on;

    % (E-2) 전압 Error
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_cc(idx), 'k', 'LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [F] CC SOC 기반 모델 : 100초 간격 줌인 (2×1 Subplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval100 = 100;
numWindows100 = floor(maxTime / interval100);

for iWin = 1:numWindows100
    tStart = (iWin-1)*interval100;
    tEnd   = iWin*interval100;
    idx = (time_s >= tStart) & (time_s < tEnd);

    figure('Name',...
        ['[CC, 100s 구간 #' num2str(iWin) '] ' ...
         '(' num2str(tStart) 's ~ ' num2str(tEnd) 's)'], ...
         'NumberTitle','off');
    c_local = lines(3);

    % (F-1) 전압 & 전류
    subplot(2,1,1);
    yyaxis left
    plot(time_s(idx), cellVoltage_meas(idx), 'Color', c_local(2,:), 'LineWidth',1.2); hold on;
    plot(time_s(idx), V_est_cc(idx),        'Color', c_local(3,:), 'LineWidth',1.2);
    ylabel('Cell Voltage (V)');
    ax = gca; ax.YColor = 'r';

    yyaxis right
    plot(time_s(idx), cellCurrent(idx), 'Color', c_local(1,:), 'LineWidth',1.2);
    ylabel('Cell Current (A)');
    ax = gca; ax.YColor = 'b';

    legend('V_{data}','V_{model}','I_{data}','Location','best');
    xlabel('Time (s)');
    grid on;

    % (F-2) 전압 Error
    subplot(2,1,2);
    plot(time_s(idx), cellVoltage_meas(idx) - V_est_cc(idx), 'k', 'LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('Voltage Error (V)');
    grid on;
end

%% (Z) BMS SOC vs CC SOC : 한 그림에 비교 플롯
figure('Name','[SOC 비교: BMS vs CC]','NumberTitle','off');
c_local = lines(3);
plot(time_s, SOC_bms, 'Color', c_local(1,:), 'LineWidth', 1.5); hold on;
plot(time_s, SOC_cc,  'Color', c_local(2,:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('SOC (%)');
legend('BMS SOC','CC SOC','Location','best');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (아래는 로컬 함수 정의)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = computeRMSE_1RC_tau(X, SOC_t, I_cell, dt, V_meas, socOCV, ocvCellVoltage)
    V_est = modelVoltage_1RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage);
    cost  = sqrt(mean((V_meas - V_est).^2));
end

function V_est = modelVoltage_1RC_tau(X, SOC_t, I_cell, dt, socOCV, ocvCellVoltage)
    % X = [R0, R1, tau1]
    R0   = X(1);
    R1   = X(2);
    tau1 = X(3);

    N = length(SOC_t);
    V_est = zeros(N,1);
    % (주의) Vrc1=0; 선언을 제거하고, 첫 샘플에서 직접 계산

    for k = 1:N
        % (1) OCV (SOC -> 전압 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear','extrap');

        % (2) R0 전압강하
        IR0 = R0 * I_cell(k);

        % (3) RC 항 계산
        alpha1 = exp(-dt(k)/tau1);
        if k == 1
            % 첫 샘플(k=1): 초기 RC 전압 설정
            Vrc1 = R1 * I_cell(k) * (1 - alpha1);
        else
            % 이후(k>1)
            Vrc1 = Vrc1*alpha1 + R1*(1 - alpha1)*I_cell(k);
        end

        % (4) 최종 모델 전압
        V_est(k) = OCV_now - IR0 - Vrc1;
    end
end

