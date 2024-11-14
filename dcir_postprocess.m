% 파일 경로 가져오기

%% 1. file load


%% 2. chg, dis step


%% 3. soc 계산

%cumtrapz


%% 4. soc = 0, 0.5, 1인 부분에 대해서 전압,전류 그래프 그리기

clc; clear; close all;

data_folder = 'C:\Users\deu04\OneDrive\문서\MATLAB\dcir';
save_path = data_folder;
I_1C = 0.00382; %[A]
id_cfa = 2; % 1 for cathode, 2 for fullcell , 3 for anode 

% MAT 파일 가져오기
slash = filesep;
files = dir([data_folder slash '*.mat']);

% 선택할 파일의 인덱스
selected_file_index = 3; % 첫 번째 파일 선택

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

for i = 1: length(data)
    if id_cfa == 1 || id_cfa == 2
        data(i).SOC = data(i).cumsumQ/total_QC;
    elseif id_cfa == 3
        data(i).SOC = 1 - data(i).cumsumQ/total_QD;
    end
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
        t_step = data(i).t - data(i).t(1); % Subtract the initial time to start from 0
        plot(t_step, V_step, 'b');
        hold on;
    end
end
xlabel('Time (s)');
ylabel('Voltage (V)');
title(ax3, 'Voltage at SOC=0');
xlim([0 30])


ax4 = subplot(1,3,2);
for i = 1:length(data)
    if abs(data(i).SOC - 0.5) < 10^(-1) % soc=0.5
        I_step = data(i).I;
        charge_mask = (I_step < 0);
        V_step = data(i).V;
        t_step = data(i).t - data(i).t(1); % Subtract the initial time to start from 0
        plot(t_step(charge_mask), V_step(charge_mask), 'b');
        hold on;
    end
end
xlabel('Time (s)');
ylabel('Voltage (V)');
title(ax4, 'Charge Pulse at SOC=0.5');
xlim([0 30])

ax5 = subplot(1,3,3);
for i = 1:length(data)
    if abs(data(i).SOC - 1) < 1e-6 % soc=1
        V_step = data(i).V;
        t_step = data(i).t - data(i).t(1); % Subtract the initial time to start from 0
        plot(t_step, V_step, 'b');
        hold on;
    end
end
xlabel('Time (s)');
ylabel('Voltage (V)');
title(ax5, 'Voltage at SOC=1');
xlim([0 30])
