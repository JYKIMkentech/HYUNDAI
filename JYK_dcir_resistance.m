% 파일 경로 가져오기
clc; clear; close all;

data_folder = 'C:\Users\deu04\OneDrive\문서\MATLAB\dcir';
save_path = data_folder;
I_1C = 0.00382; %[A]

% MAT 파일 가져오기
slash = filesep;
files = dir([data_folder slash '*.mat']);

for i = 1:length(files)
   fullpath_now = [data_folder slash files(i).name];% path for i-th file in the folder
   load(fullpath_now);
   data(1)= [];

end
% 충전, 방전 스텝(필드) 구하기 

step_chg = [];
step_dis = [];

for i = 1:length(data)
    % type 필드가 C인지 확인
    if strcmp(data(i).type, 'C')
        % C가 맞으면 idx 1 추가
        step_chg(end+1) = i;
    % type 필드가 D인지 확인
    elseif strcmp(data(i).type, 'D')

        % 맞으면 idx 1 추가
        step_dis(end+1) = i;
    end
end



% STEP 내부에서의 전하량 구하기

for j = 1:length(data)
     %calculate capacities
     data(j).Q = trapz(data(j).t,data(j).I)/3600; %[Ah]
     data(j).cumQ = cumtrapz(data(j).t,data(j).I)/3600; %[Ah]
     

     % data(j).cumQ = abs(cumtrapz(data(j).t,data(j).I))/3600; %[Ah]
     
end

% Total QC, QD값 구하기 ( 전체 전하량 구하기) 
total_QC = sum(abs([data(step_chg).Q]));  % charge 상태 전체 Q값
total_QD = sum(abs([data(step_dis).Q])); % discharge 상태 전체 Q값



% cumsumQ 필드 추가
for i = 1:length(data)
    if i == 1
        data(i).cumsumQ = data(i).cumQ;
    else
        data(i).cumsumQ = data(i-1).cumsumQ(end) + data(i).cumQ;
    end
end

for i = 1 : length(data)
    data(i).SOC = data(i).cumsumQ/total_QC;
end


% 큰 I 가지는 index 추출
BigI = [];

for i = 1:length(data)
    if abs(data(i).I) > (1/3 * I_1C)
       BigI = [BigI , i];
    end
end

% I의 평균을 필드에 저장하기 

for i = 1:length(data)
    data(i).avgI = mean(data(i).I);
end

% V 변화량 구하기
for i = 1 : length(data)
    if i == 1
       data(i).deltaV = zeros(size(data(i).V));
    else
       data(i).deltaV = data(i).V() - data(i-1).V(end);
    end
end

% Resistance 구하기 
for i = 1 : length(data)
    if data(i).avgI == 0
        data(i).R = zeros(size(data(i).V));
    else 
        data(i).R = (data(i).deltaV / data(i).avgI) .* ones(size(data(i).V));
    end
end
% plot

% BigI 부분에서만 plot 하게 만들기

% R 부분은 저항 0으로 하기, BigI에 해당하는 data번호만 그리기

plot(data(6).t, data(6).R)

xlabel('time (sec)')
ylabel('Resistance')

% x 값이 30초일 때의 y 값을 얻기
x_value = 7201;
y_value = interp1(data(2).t, data(2).R, x_value);

disp(y_value); % 결과 출력

%0.01 sec 에서 Resistance 
for i = 1:length(BigI)
    x_001 = data(BigI(i)).t(1) + 0.01;
    data(BigI(i)).R001s = interp1(data(BigI(i)).t, data(BigI(i)).R , x_001);
end

% 1s , 10s, 30s 에서 Resistance 
for i = 1:length(BigI)
   data(BigI(i)).R1s = data(BigI(i)).R(11);
   data(BigI(i)).R10s = data(BigI(i)).R(56);
   data(BigI(i)).R30s = data(BigI(i)).R(end);
end

% BigI에서 charge 상태까지의 step 얻기
BigIC = BigI(BigI < step_chg(end));
BigID = BigI(BigI >= step_chg(end));


% 10s

% 30 s

% 데이터의 차이는 0.1초 = 100ms

% SOC-Resistance 그래프 그리기

% 각각의 Resistance에 대응되는 시간 - SOC 지정하기 





% 0.001s

SOC001s = [];
R001s = [];
for i = 1:length(BigIC)
    SOC001s = [SOC001s, data(BigIC(i)).SOC(2)];
    R001s = [R001s, data(BigIC(i)).R001s];
end


% 1s

SOC1s = [];
R1s = [];
for i = 1:length(BigIC)
    SOC1s = [SOC1s, data(BigIC(i)).SOC(11)];
    R1s = [R1s, data(BigIC(i)).R(11)];
end


% 10s

SOC10s = [];
R10s = [];
for i = 1:length(BigIC)
    SOC10s = [SOC10s, data(BigIC(i)).SOC(56)];
    R10s = [R10s, data(BigIC(i)).R10s];
end



% 30s

SOC30s = [];
R30s = [];
for i = 1:length(BigIC)
    SOC30s = [SOC30s, data(BigIC(i)).SOC(end)];
    R30s = [R30s, data(BigIC(i)).R(end)];
end

% spline을 사용하여 점들을 이어주기
smoothed_SOC_001s = linspace(min(SOC001s), max(SOC001s), 100);
smoothed_R_001s = spline(SOC001s, R001s, smoothed_SOC_001s);

smoothed_SOC_1s = linspace(min(SOC1s), max(SOC1s), 100); 
smoothed_R_1s = spline(SOC1s, R1s, smoothed_SOC_1s); 

smoothed_SOC_10s = linspace(min(SOC10s), max(SOC10s), 100);
smoothed_R_10s = spline(SOC10s, R10s, smoothed_SOC_10s);

smoothed_SOC_30s = linspace(min(SOC30s), max(SOC30s), 100); 
smoothed_R_30s = spline(SOC30s, R30s, smoothed_SOC_30s); 

% 그래프 그리기
figure;
hold on;
plot(SOC001s, R001s, 'o');
plot(smoothed_SOC_001s, smoothed_R_001s);
plot(SOC1s, R1s, 'o');
plot(smoothed_SOC_1s, smoothed_R_1s);
plot(SOC10s, R10s, 'o');
plot(smoothed_SOC_10s, smoothed_R_10s);
plot(SOC30s, R30s, 'o');
plot(smoothed_SOC_30s, smoothed_R_30s);
hold off;

xlabel('SOC');
ylabel('Resistance (\ohm cm^2)', 'fontsize', 12);
title('SOC vs Resistance');
legend('100ms', '100ms (line)', '1s', '1s (line)', '10s', '10s (line)', '30s', '30s (line)'); 


