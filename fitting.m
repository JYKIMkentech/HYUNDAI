load("dcir_fit.mat");

% 모델 함수 정의
function voltage = model_func(time, R1, R2, C, a)
    voltage = I * R1 * (R1 + R2) ./ (R1 + R2 .* exp(-a .* time ./ (R1 * C)));
end


