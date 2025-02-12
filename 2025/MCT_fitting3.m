%% script_2RC_TauSurface_paramArray.m
% (목적) tau1, tau2를 격자(grid)로 고정하고, R0,R1,R2만 최적화하여 RMSE를 계산
%        -> (tau1, tau2) 에 대한 Cost(RMSE) Surface를 그림
%        -> paramArray에 초기값, 최적화값, RMSE 등을 모두 저장

clear; clc; close all;

%% 1) 사용자 입력 (MCT 번호)
mctNumber = input('피팅할 MCT 번호를 입력하세요 (1~6): ');

%% 2) 코드1에서 저장한 결과 불러오기
%    (mctCellData{mctNumber} 안에 Time_s, PackVoltage_V, Current_A 등이 있다고 가정)
load('MCT_Results.mat', ...
     'mctCellData','OCVMCT','dataOCV', ...
     'socOCV','ocvCellVoltage','uSocOCV','uCellVoltage');

%% 3) 사용할 배터리 총용량 (예: Q_batt = 56.2396 [Ah])
Q_batt = 56.2396;  

%% 4) fmincon 옵션 설정 (R0,R1,R2만 최적화)
%    --> 여기서 x0(1)=R0, x0(2)=R1, x0(3)=R2가 초기값
x0 = [0.001, 0.0005, 0.0005];  % 초기 guess [R0, R1, R2]
lb = [0, 0, 0];
ub = [inf, inf, inf];

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                       'MaxIterations',1000,'MaxFunctionEvaluations',5000);

%% 5) 해당 MCT 데이터 로드
dataMCT = mctCellData{mctNumber};

time_s      = dataMCT.Time_s;
packVoltage = dataMCT.PackVoltage_V;
packCurrent = dataMCT.Current_A;   % (+)가 방전 전류

% 192 직렬, 2 병렬 가정
cellVoltage_meas = packVoltage / 192;  % 셀 전압
cellCurrent      = packCurrent / 2;    % 셀 전류

%% 6) 초기 SOC 계산
cellVoltage_init = mean(cellVoltage_meas(1:5));  
SOC0 = interp1(ocvCellVoltage, socOCV, cellVoltage_init, 'linear', 'extrap');

% 시계열 간격
dt = [0; diff(time_s)];

% 시계열 SOC(t) 계산
charge_integral = cumtrapz(time_s, cellCurrent);  % [A·s]
SOC_t = SOC0 - (charge_integral / (Q_batt*3600))*100;  % [%]

%% 7) tau1, tau2 격자 생성
tau1_candidates = linspace(0.1, 10, 21);   % (예) 0.1 ~ 10, 총 21점
tau2_candidates = linspace(20, 80, 21);   % (예) 20 ~ 80, 총 21점

% costSurface(i,j) = RMSE at (tau1_candidates(i), tau2_candidates(j))
costSurface = zeros(length(tau1_candidates), length(tau2_candidates));

% === (중요) paramArray에 [초기 R0,R1,R2, 초기 tau1,tau2, 최적 R0,R1,R2, 최적 tau1,tau2, RMSE] 저장 ===
% 크기: (길이(tau1_candidates), 길이(tau2_candidates), 11)
paramArray = zeros(length(tau1_candidates), length(tau2_candidates), 11);

%% 8) 2중 for문으로 tau1, tau2 고정 후 R0,R1,R2 최적화
for i = 1:length(tau1_candidates)
    for j = 1:length(tau2_candidates)
        tau1_now = tau1_candidates(i);
        tau2_now = tau2_candidates(j);

        % costFunc 정의 (R0,R1,R2만 최적화)
        costFunc = @(x) computeRMSE_2RC_3param(... 
            x, tau1_now, tau2_now, ...
            SOC_t, cellCurrent, dt, ...
            cellVoltage_meas, socOCV, ocvCellVoltage);

        % fmincon 실행
        [xOpt, fVal] = fmincon(costFunc, x0, [], [], [], [], lb, ub, [], options);

        % RMSE (cost)
        costSurface(i,j) = fVal;

        % ===== paramArray에 저장 =====
        % 1) 초기 R0, 2) 초기 R1, 3) 초기 R2
        % 4) 초기 tau1(= tau1_now), 5) 초기 tau2(= tau2_now)
        % 6) 피팅된 R0(xOpt(1)), 7) 피팅된 R1(xOpt(2)), 8) 피팅된 R2(xOpt(3))
        % 9) 피팅된 tau1(= tau1_now 고정), 10) 피팅된 tau2(= tau2_now 고정)
        % 11) RMSE
        paramArray(i,j,:) = [ ...
            x0(1), x0(2), x0(3), ...
            tau1_now, tau2_now, ...
            xOpt(1), xOpt(2), xOpt(3), ...
            tau1_now, tau2_now, ...
            fVal ...
        ];
    end
end

%% 9) costSurface에서 최소값 찾기
[bestRMSE, idxMin] = min(costSurface(:));
[best_i, best_j]   = ind2sub(size(costSurface), idxMin);

best_tau1  = tau1_candidates(best_i);
best_tau2  = tau2_candidates(best_j);

% paramArray에서 최적화된 R0,R1,R2 뽑아오기
best_R0 = paramArray(best_i, best_j, 6); 
best_R1 = paramArray(best_i, best_j, 7);
best_R2 = paramArray(best_i, best_j, 8);

fprintf('\n=== Best RMSE = %.4f ===\n', bestRMSE);
fprintf('   tau1* = %.3f, tau2* = %.3f\n', best_tau1, best_tau2);
fprintf('   [R0,R1,R2] = [%.4f, %.4f, %.4f]\n', best_R0, best_R1, best_R2);

%% [추가] 사용자 정의 "초기 추정" (tau1_init, tau2_init)에 대한 결과 확인
tau1_init = 6;   % 예: 6초
tau2_init = 30;  % 예: 30초

% (1) 격자 중에서 init값과 가장 가까운 인덱스 찾기
[~, idx_i_guess] = min(abs(tau1_candidates - tau1_init));
[~, idx_j_guess] = min(abs(tau2_candidates - tau2_init));

% (2) 해당 인덱스에서의 RMSE 및 (R0,R1,R2)
costAtGuess  = costSurface(idx_i_guess, idx_j_guess);

% paramArray에서 R0,R1,R2 가져오기
paramAtGuess = squeeze(paramArray(idx_i_guess, idx_j_guess, :));  
% -> [initR0, initR1, initR2, initTau1, initTau2, optR0, optR1, optR2, optTau1, optTau2, RMSE]

fprintf('\n--- [Initial guess] tau1=%.2f, tau2=%.2f ---\n', tau1_init, tau2_init);
fprintf('    -> RMSE = %.4f\n', costAtGuess);
fprintf('    -> [R0,R1,R2] = [%.4f, %.4f, %.4f]\n', ...
    paramAtGuess(6), paramAtGuess(7), paramAtGuess(8));

%% 10) Surf 그래프 (3D) 그리기 + 최솟값 지점, 초기 추정 지점 표시
[T1, T2] = meshgrid(tau2_candidates, tau1_candidates);  
% ※ surf() 시 X->열방향(tau2), Y->행방향(tau1) 순서이므로 주의

figure('Name','RMSE Surface','NumberTitle','off');
surfH = surf(T1, T2, costSurface, ...
    'EdgeColor','none','FaceColor','interp', ...  % 보기 좋게
    'DisplayName','RMSE Surface');                
shading interp; colorbar; grid on;
xlabel('\tau_2'); ylabel('\tau_1'); zlabel('RMSE (V)');
title(['MCT-' num2str(mctNumber) ' : RMSE Surface']);
hold on;

% (1) 최소값 지점(빨간 동그라미)
minH = plot3(best_tau2, best_tau1, bestRMSE, ...
    'ro','MarkerSize',10,'LineWidth',1.5, ...
    'DisplayName', ...
    sprintf(['Min RMSE point \n(Opt \\tau_1=%.2f, \\tau_2=%.2f)'], ...
             best_tau1, best_tau2));

% (2) 초기 추정값 지점(파란 동그라미)
guessH = plot3(tau2_init, tau1_init, costAtGuess, ...
    'bo','MarkerSize',10,'LineWidth',1.5, ...
    'DisplayName', ...
    sprintf(['Initial guess \n(\\tau_1=%.2f, \\tau_2=%.2f)'], ...
             tau1_init, tau2_init));

legend([surfH, minH, guessH], 'Location','best');

%% (옵션) Contour(2D)로도 표시
figure('Name','RMSE Contour (tau1 vs tau2)','NumberTitle','off');
contourf(tau2_candidates, tau1_candidates, costSurface, 20, 'LineColor','none');
hold on;

% (1) 최소값 지점 표시 (2D 평면 빨간 점)
plot(best_tau2, best_tau1, 'ro','MarkerSize',10,'LineWidth',1.5);

% (2) 초기값 지점 표시 (2D 평면 파란 점)
plot(tau2_init, tau1_init, 'bo','MarkerSize',10,'LineWidth',1.5);

xlabel('\tau_2'); ylabel('\tau_1');
title(['MCT-' num2str(mctNumber) ' : RMSE Contour']);
colorbar; grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (아래는 로컬 함수 정의부)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = computeRMSE_2RC_3param(x, tau1, tau2, ...
                                       SOC_t, I_cell, dt, ...
                                       V_meas, socOCV, ocvCellVoltage)
    % x = [R0, R1, R2]
    R0 = x(1);
    R1 = x(2);
    R2 = x(3);

    % 모델 전압
    V_est = modelVoltage_2RC_3param(R0, R1, R2, tau1, tau2, ...
                                    SOC_t, I_cell, dt, ...
                                    socOCV, ocvCellVoltage);

    % RMSE
    cost = sqrt(mean((V_meas - V_est).^2));
end

function V_est = modelVoltage_2RC_3param(R0, R1, R2, tau1, tau2, ...
                                         SOC_t, I_cell, dt, ...
                                         socOCV, ocvCellVoltage)
    % 2개의 RC (전압 항) 초기값
    Vrc1 = 0;
    Vrc2 = 0;

    N = length(SOC_t);
    V_est = zeros(N,1);

    for k = 1:N
        % (1) OCV (SOC 기반 보간)
        OCV_now = interp1(socOCV, ocvCellVoltage, SOC_t(k), 'linear', 'extrap');

        % (2) R0 전압 강하
        IR0 = R0 * I_cell(k);

        % (3) RC 업데이트
        if k > 1
            alpha1 = exp(-dt(k)/tau1);
            alpha2 = exp(-dt(k)/tau2);

            Vrc1 = Vrc1*alpha1 + R1*(1 - alpha1)*I_cell(k);
            Vrc2 = Vrc2*alpha2 + R2*(1 - alpha2)*I_cell(k);
        end

        % (4) 최종 전압
        V_est(k) = OCV_now - IR0 - Vrc1 - Vrc2;
    end
end

