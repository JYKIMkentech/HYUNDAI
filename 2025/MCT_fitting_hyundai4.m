clc; clear; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 코드1에서 저장한 결과 불러오기
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

%% 3) fmincon 옵션 설정 (R0, R1만 최적화)
%    x = [R0, R1]
x0 = [0.001, 0.0005];  % 초기값
lb = [0, 0];
ub = [inf, inf];

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxIterations',1000,'MaxFunctionEvaluations',5000);

%% 4) 해당 MCT 데이터 로드
dataMCT       = mctCellData{mctNumber};
time_s        = dataMCT.Time_s;          % (s)
packVoltage   = dataMCT.PackVoltage_V;   % (V)
packCurrent   = dataMCT.Current_A;       % (양수=방전)
socDecimal    = dataMCT.SOC_decimal;     % (소수점 SOC)
socInteger    = dataMCT.SOC_integer;     % (정수 SOC)

% (A) BMS SOC
SOC_bms = socInteger + socDecimal;  % [%]

% (B) Cell 전압/전류 변환 (예: 192직렬, 2병렬)
cellVoltage_meas = packVoltage / 192;  
cellCurrent      = packCurrent / 2;    

% (C) 시간간격 dt
dt = [0; diff(time_s)];

%% 5) tau1 후보군 생성
tau1_candidates = linspace(0, 500, 501);  % 0~200 사이를 201점으로 균등분할
costArray       = zeros(length(tau1_candidates), 1);

% (중요) 각 tau1에서의 결과 저장: [초기R0, 초기R1, 고정tau1, 최적R0, 최적R1, 고정tau1, RMSE]
paramArray      = zeros(length(tau1_candidates), 7);

%% 6) for문: tau1을 고정하고 R0,R1만 최적화
for i = 1:length(tau1_candidates)
    tau1_now = tau1_candidates(i);

    % costFunc 정의
    costFunc = @(x) computeRMSE_1RC_2param( ...
        x, tau1_now, ...
        SOC_bms, cellCurrent, dt, ...
        cellVoltage_meas, socOCV, ocvCellVoltage);

    % fmincon 실행
    [xOpt, fVal] = fmincon(costFunc, x0, [], [], [], [], lb, ub, [], options);

    % RMSE 저장
    costArray(i) = fVal;

    % 결과 기록
    % [초기 R0, 초기 R1, 고정 tau1, 최적 R0, 최적 R1, 고정 tau1, RMSE]
    paramArray(i,:) = [x0(1), x0(2), tau1_now, xOpt(1), xOpt(2), tau1_now, fVal];
end

%% 7) costArray에서 최소값 찾기
[bestRMSE, idxMin] = min(costArray);
best_tau1 = tau1_candidates(idxMin);
best_R0   = paramArray(idxMin, 4);
best_R1   = paramArray(idxMin, 5);

fprintf('\n=== Best RMSE = %.4f ===\n', bestRMSE);
fprintf('   tau1* = %.3f\n', best_tau1);
fprintf('   [R0,R1] = [%.4f, %.4f]\n', best_R0, best_R1);

%% 8) RMSE vs tau1 플롯
figure('Name','RMSE vs tau1','NumberTitle','off');
plot(tau1_candidates, costArray, 'b-o','LineWidth',1.5,'MarkerSize',6); hold on;
plot(best_tau1, bestRMSE, 'ro','MarkerSize',10,'LineWidth',1.5);
xlabel('\tau_1'); ylabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : RMSE vs \tau_1']);
grid on;
legend('RMSE','Minimum RMSE','Location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (아래는 로컬 함수 정의)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = computeRMSE_1RC_2param(x, tau1, ...
                                       SOC_bms, I_cell, dt, ...
                                       V_meas, socOCV, ocvCellVoltage)
    % x = [R0, R1]
    R0 = x(1);
    R1 = x(2);

    % 모델 전압
    V_est = modelVoltage_1RC_2param(R0, R1, tau1, ...
                                    SOC_bms, I_cell, dt, ...
                                    socOCV, ocvCellVoltage);

    % RMSE
    cost = sqrt(mean((V_meas - V_est).^2));
end

function V_est = modelVoltage_1RC_2param(R0, R1, tau1, ...
                                         SOC_bms, I_cell, dt, ...
                                         socOCV, ocvCellVoltage)
    N   = length(SOC_bms);
    Vrc = 0;                % RC 항 초기값
    V_est = zeros(N, 1);    % 전체 시점에 대한 모델 전압

    for k = 1:N
        % (1) OCV (SOC 기반 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_bms(k), 'linear', 'extrap');

        % (2) R0 전압 강하
        IR0 = R0 * I_cell(k);

        % (3) RC 전압 업데이트
        if k > 1
            alpha = exp(-dt(k)/tau1);
            Vrc   = Vrc*alpha + R1*(1 - alpha)*I_cell(k);
        end

        % (4) 최종 전압
        V_est(k) = OCV_now - IR0 - Vrc;
    end
end

