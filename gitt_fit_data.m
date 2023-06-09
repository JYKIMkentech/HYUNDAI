function fit_data()
    % 데이터 로드
    load('gitt_fit.mat');
    deltaV_exp = data(22).deltaV;
    time_exp = data(22).t;
    
    % 최적화를 위한 초기 추정값
    R1 = 12.5032;
    R2 = 24.4340;
    C = 0.00806;
    A = data(22).V(56) - data(22).V(133)
    B = 1200

    initial_guess = [R1, R2, C, A, B];
    
    % fmincon을 사용하여 최적화 수행
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 100);
    [opt_params, rms] = fmincon(@(params) cost_function(params, time_exp, deltaV_exp), ...
        initial_guess, [], [], [], [], [0, 0, 0, 0, 0], [], [], options); % Added two zeros to lower bound constraint to match the size of initial_guess
    
    % 최적화된 파라미터 출력
    disp("Optimized Parameters:");
    disp("R1: " + opt_params(1));
    disp("R2: " + opt_params(2));
    disp("C: " + opt_params(3));
    disp("A: " + opt_params(4));
    disp("B: " + opt_params(5));
    
    % 최적화된 파라미터를 사용하여 모델 예측
    voltage_model = model_func(time_exp, opt_params(1), opt_params(2), opt_params(3), opt_params(4), opt_params(5));
    
    % 데이터와 모델 결과를 그래프로 플롯
    plot(time_exp, deltaV_exp, 'b-', time_exp, voltage_model, 'r--');
    legend('실험 데이터', '모델 결과');
    xlabel('시간');
    ylabel('전압');
    title('실험 데이터와 모델 결과');
end


function cost = cost_function(params, time, deltaV)
    R1 = params(1);
    R2 = params(2);
    C = params(3);
    A = params(4);
    B = params(5);
    
    % 모델 함수를 사용하여 예측 전압 계산
    voltage_model = model_func(time, R1, R2, C, A, B);
    
    % RMS 오차 계산
    error = deltaV - voltage_model;
    cost = sqrt(mean(error.^2));

end

% 모델 함수 정의
function voltage = model_func(time, R1, R2, C, A, B)
    I = 0.0038;
    voltage = zeros(size(time)); % 전압을 저장할 벡터 초기화
    for i = 1:length(time)
        t = time(i); % 각 시간 단계에 대해 전압 계산
        voltage(i) = I * R1 * (R1 + R2 + A * (1-sqrt(t/B))) / (R1 + (R2 + A * (1-sqrt(t/B))) * exp((-R1/(R2+A * (1-sqrt(t/B))) + 1) * t / (R1 * C)));
    end
end
