function fit_data()
    % 데이터 로드
    load('dcir_fit.mat');
    deltaV_exp = data(22).deltaV;
    time_exp = data(22).t;
    
    % 최적화를 위한 초기 추정값
    initial_guess = [12.5032, 24.4340, 0.00806];
    
    % fmincon을 사용하여 최적화 수행
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 100);
    [opt_params, rms] = fmincon(@(params) cost_function(params, time_exp, deltaV_exp), ...
        initial_guess, [], [], [], [], [0, 0, 0], [], [], options);
    
    % 최적화된 파라미터 출력
    disp("Optimized Parameters:");
    disp("R1: " + opt_params(1));
    disp("R2: " + opt_params(2));
    disp("C: " + opt_params(3));
    
    % 최적화된 파라미터를 사용하여 모델 예측
    voltage_model = model_func(time_exp, opt_params(1), opt_params(2), opt_params(3));
    
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
    
    % 모델 함수를 사용하여 예측 전압 계산
    voltage_model = model_func(time, R1, R2, C);
    
    % RMS 오차 계산
    error = deltaV - voltage_model;
    cost = sqrt(mean(error.^2));
end

function voltage = model_func(time, R1, R2, C)
    I = 0.0038;  % 초기에 정의한 I 값을 사용
    a = (R2/R1)+1;  % 초기에 정의한 a 값을 사용
    
    voltage = I * R1 * (R1 + R2) ./ (R1 + R2 .* exp(-a .* time ./ (R1 * C)));
end