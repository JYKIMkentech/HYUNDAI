% 파일 경로 가져오기
clc; clear; close all;

data_folder = 'C:\Users\deu04\OneDrive\문서\MATLAB\dcir';
save_path = data_folder;
I_1C = 0.00382; %[A]
id_cfa = 1; % 1 for cathode, 2 for fullcell , 3 for anode 

% MAT 파일 가져오기
slash = filesep;
files = dir([data_folder slash '*.mat']);

% 선택할 파일의 인덱스
selected_file_index = 2; % 첫 번째 파일 선택

% 선택한 파일 load
fullpath_now = [data_folder slash files(selected_file_index).name];
load(fullpath_now);
data(1) = [];

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
     data(j).Q = abs(trapz(data(j).t,data(j).I))/3600; %[Ah]
     data(j).cumQ = abs(cumtrapz(data(j).t,data(j).I))/3600; %[Ah]
     

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
    if id_cfa == 1 || id_cfa == 2
        if id_cfa == 1
            data(i).SOC = data(i).cumsumQ/total_QC; % Cathode
        elseif id_cfa == 2
            data(i).SOC = data(i).cumsumQ/total_QC; % FCC
        end
        % 큰 I 가지는 index 추출
        BigI = [];
        for i = 1:length(data)
            if abs(data(i).I) > (1/3 * I_1C)
               BigI = [BigI , i];
            end
        end
        
        if id_cfa == 1 || id_cfa == 2
            % BigIC, BigID 계산
            BigIC = BigI(BigI < step_chg(end));
            BigID = BigI(BigI >= step_chg(end));
        end
    elseif id_cfa == 3
        BigI = [];
        for i = 1:length(data)
            data(i).SOC = 1 - data(i).cumsumQ/total_QD;
            if abs(data(i).I) > (1/3 * I_1C)
               BigI = [BigI , i];
               
            end
        end
        % BigI 계산
        BigI = BigI;
       
    else
        error('Invalid id_cfa value. Please choose 1 for cathode, 2 for FCC, or 3 for anode.');
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
% -- CATHOD, FCC 에서는 BIGIC,D 구간 얻기 -- %
% BigIC = BigI(BigI < step_chg(end));
% BigID = BigI(BigI >= step_dis(end));


% 10s
% 30 s
% 데이터의 차이는 0.1초 = 100ms
% SOC-Resistance 그래프 그리기
% 각각의 Resistance에 대응되는 시간 - SOC 지정하기

% 0.001s
% CATHODE, FCC = BIGIC 데이터 확인
% ANDOE = BIGI 데이터 확인
SOC001s = [];
R001s = [];
SOC1s = [];
R1s = [];
SOC10s = [];
R10s = [];
SOC30s = [];
R30s = [];

if id_cfa == 1 || id_cfa == 2
    for i = 1:length(BigIC)
        SOC001s = [SOC001s, data(BigIC(i)).SOC(2)];
        R001s = [R001s, data(BigIC(i)).R001s];
        SOC1s = [SOC1s, data(BigIC(i)).SOC(11)];
        R1s = [R1s, data(BigIC(i)).R1s];
        SOC10s = [SOC10s, data(BigIC(i)).SOC(56)];
        R10s = [R10s, data(BigIC(i)).R10s];
        SOC30s = [SOC30s, data(BigIC(i)).SOC(end)];
        R30s = [R30s, data(BigIC(i)).R(end)];
    end
elseif id_cfa == 3
    for i = 1:length(BigI)
        SOC001s = [SOC001s, data(BigI(i)).SOC(2)];
        R001s = [R001s, data(BigI(i)).R001s];
        SOC1s = [SOC1s, data(BigI(i)).SOC(11)];
        R1s = [R1s, data(BigI(i)).R1s];
        SOC10s = [SOC10s, data(BigI(i)).SOC(56)];
        R10s = [R10s, data(BigI(i)).R10s];
        SOC30s = [SOC30s, data(BigI(i)).SOC(end)];
        R30s = [R30s, data(BigI(i)).R(end)];
    end
end

% % 1s
% SOC1s = [];
% R1s = [];
% for i = 1:length(BigIC)
%     if id_cfa == 1 || id_cfa == 2
%         SOC1s = [SOC1s, data(BigIC(i)).SOC(11)];
%         R1s = [R1s, data(BigIC(i)).R1s];
%     elseif id_cfa == 3
%         SOC1s = [SOC1s, data(BigI(i)).SOC(11)];
%         R1s = [R1s, data(BigI(i)).R1s];
%     end
% end
% 
% % 10s
% SOC10s = [];
% R10s = [];
% for i = 1:length(BigIC)
%     if id_cfa == 1 || id_cfa == 2
%         SOC10s = [SOC10s, data(BigIC(i)).SOC(56)];
%         R10s = [R10s, data(BigIC(i)).R10s];
%     elseif id_cfa == 3
%         SOC10s = [SOC10s, data(BigI(i)).SOC(56)];
%         R10s = [R10s, data(BigI(i)).R10s];
%     end
% end

% % 30s
% SOC30s = [];
% R30s = [];
% for i = 1:length(BigIC)
%     if id_cfa == 1 || id_cfa == 2
%         SOC30s = [SOC30s, data(BigIC(i)).SOC(end)];
%         R30s = [R30s, data(BigIC(i)).R(end)];
%     elseif id_cfa == 3
%         SOC30s = [SOC30s, data(BigI(i)).SOC(end)];
%         R30s = [R30s, data(BigI(i)).R(end)];
%     end
% end
% spline을 사용하여 점들을 부드럽게 이어주기
% Smoothed SOC-Resistance curve for SOC001s
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
ylabel('Resistance (\Omega )', 'fontsize', 12);
title('SOC vs Resistance');
legend('100ms', '100ms (line)', '1s', '1s (line)', '10s', '10s (line)', '30s', '30s (line)' ,'Location', 'northwest'); 
xlim([0 1])


