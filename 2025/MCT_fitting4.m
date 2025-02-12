%% script_1RC_TauSurface.m
% (목적) tau1를 격자(grid)로 고정하고, R0, R1만 최적화하여 RMSE를 계산
%        -> (tau1)에 대한 Cost(RMSE) 곡선을 그림

clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 코드1에서 저장한 결과 불러오기
%    (mctCellData{mctNumber} 안에 Time_s, PackVoltage_V, Current_A 등의 필드가 있다고 가정)
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

%% 3) 사용할 배터리 총용량 (예: Q_batt = 56.2396 [Ah])
Q_batt = 56.2396;  

%% 4) fmincon 옵션 설정 (R0,R1만 최적화)
x0 = [0.001, 0.0005];  % 초기값 [R0, R1]
lb = [0, 0];
ub = [inf, inf];

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                       'MaxIterations',1000,'MaxFunctionEvaluations',5000);

%% 5) 해당 MCT 데이터 로드
dataMCT = mctCellData{mctNumber};

time_s      = dataMCT.Time_s;
packVoltage = dataMCT.PackVoltage_V;
packCurrent = dataMCT.Current_A;   % (+)가 방전 전류

% 예) 192 직렬, 2 병렬
cellVoltage_meas = packVoltage / 192;  % 셀 전압
cellCurrent      = packCurrent / 2;    % 셀 전류

%% 6) 초기 SOC 계산
% 첫 몇 개 샘플 평균 전압 -> 이를 OCV 테이블 기반 보간
cellVoltage_init = mean(cellVoltage_meas(1:5));  
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear', 'extrap');

% 시계열 간격
dt = [0; diff(time_s)];

% 시계열 SOC(t) 계산
charge_integral = cumtrapz(time_s, cellCurrent);  % [A·s]
SOC_t = SOC0 - (charge_integral / (Q_batt*3600))*100;  % [%]

%% 7) tau1 후보군 생성 (0부터 20까지 총 21개)
tau1_candidates = linspace(0, 50, 51);

% costArray(i) = RMSE at tau1_candidates(i)
costArray = zeros(length(tau1_candidates), 1);

% (옵션) 각 tau1에 대해 최적화된 [R0,R1]도 저장
paramArray = zeros(length(tau1_candidates), 2);

%% 8) for문으로 tau1 고정 후 R0,R1 최적화
for i = 1:length(tau1_candidates)
    tau1_now = tau1_candidates(i);

    % costFunc 정의 (R0,R1만 최적화)
    costFunc = @(x) computeRMSE_1RC_2param(...
        x, tau1_now, ...
        SOC_t, cellCurrent, dt, ...
        cellVoltage_meas, socOCV, ocvCellVoltage);

    % fmincon 실행
    [xOpt, fVal] = fmincon(costFunc, x0, [], [], [], [], lb, ub, [], options);

    % RMSE (cost)
    costArray(i) = fVal;

    % 최적화된 파라미터 저장
    paramArray(i,:) = xOpt;
end

%% 9) costArray에서 최소값 찾기
[bestRMSE, idxMin] = min(costArray);
best_tau1 = tau1_candidates(idxMin);
best_R0R1 = paramArray(idxMin, :);  % [R0, R1]

fprintf('\n=== Best RMSE = %.4f ===\n', bestRMSE);
fprintf('   tau1* = %.3f\n', best_tau1);
fprintf('   [R0,R1] = [%.4f, %.4f]\n', best_R0R1(1), best_R0R1(2));

%% 10) RMSE vs tau1 플롯
figure('Name','RMSE vs. tau1','NumberTitle','off');
plot(tau1_candidates, costArray, 'b-o','LineWidth',1.5,'MarkerSize',6);
hold on;

% 최소값 지점 표시
plot(best_tau1, bestRMSE, 'ro','MarkerSize',10,'LineWidth',1.5);

xlabel('\tau_1');
ylabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : RMSE vs \tau_1']);
grid on;

legend('RMSE','Minimum RMSE','Location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (아래는 로컬 함수 정의부)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = computeRMSE_1RC_2param(x, tau1, ...
                                       SOC_t, I_cell, dt, ...
                                       V_meas, socOCV, ocvCellVoltage)
    % x = [R0, R1]
    R0 = x(1);
    R1 = x(2);

    % 모델 전압
    V_est = modelVoltage_1RC_2param(R0, R1, tau1, ...
                                    SOC_t, I_cell, dt, ...
                                    socOCV, ocvCellVoltage);

    % RMSE
    cost = sqrt(mean((V_meas - V_est).^2));
end

function V_est = modelVoltage_1RC_2param(R0, R1, tau1, ...
                                         SOC_t, I_cell, dt, ...
                                         socOCV, ocvCellVoltage)
    % 1개의 RC (전압 항) 초기값
    Vrc = 0;

    N = length(SOC_t);
    V_est = zeros(N,1);

    for k = 1:N
        % (1) OCV (SOC 기반 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');

        % (2) R0 전압 강하
        IR0 = R0 * I_cell(k);

        % (3) RC 업데이트
        if k > 1
            alpha = exp(-dt(k)/tau1);
            Vrc   = Vrc*alpha + R1*(1-alpha)*I_cell(k);
        end

        % (4) 최종 전압
        V_est(k) = OCV_now - IR0 - Vrc;
    end
end
