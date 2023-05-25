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




% Plot
% 전체적인 soc에 따른 전압, 전류 그래프
figure
ax1 = subplot(1,2,1);
for i = 1:length(data)
    soc_step = data(i).SOC;
    V_step = data(i).V;
    plot(soc_step, V_step, 'b');
    hold on;
end
xlabel('State of Charge (SOC)');
ylabel('Voltage (V)');
title(ax1, 'Voltage vs SOC');

ax2 = subplot(1,2,2);
for i = 1:length(data)
    soc_step = data(i).SOC;
    I_step = data(i).I;
    plot(soc_step, I_step, 'r');
    hold on;
end
xlabel('State of Charge (SOC)');
ylabel('Current (A)');
title(ax2, 'Current vs SOC');

% soc=0, 0.5, 1 일때의 전압 그래프
figure
ax3 = subplot(1,3,1);
for i = 1:length(data)
    if abs(data(i).SOC) < 1e-6 % soc=0
        V_step = data(i).V;
        plot(data(i).t, V_step, 'b');
        hold on;
    end
end
xlabel('Time (s)');
ylabel('Voltage (V)');
title(ax3, 'Voltage at SOC=0');

ax4 = subplot(1,3,2);
for i = 1:length(data)
    if abs(data(i).SOC - 0.5) < 1e-1 % soc=0.5
        I_step = data(i).I;
        charge_mask = (I_step > 0);
        V_step = data(i).V;
        plot(data(i).t(charge_mask), V_step(charge_mask), 'b');
        hold on;
    end
end
xlabel('Time (s)');
ylabel('Voltage (V)');
title(ax4, 'Charge Pulse at SOC=0.5');

ax5 = subplot(1,3,3);
for i = 1:length(data)
    if abs(data(i).SOC - 1) < 1e-6 % soc=1
        V_step = data(i).V;
        plot(data(i).t, V_step, 'b');
        hold on;
    end
end
xlabel('Time (s)');
ylabel('Voltage (V)');
title(ax5, 'Voltage at SOC=1');


