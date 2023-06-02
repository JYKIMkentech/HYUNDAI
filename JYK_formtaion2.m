clc; clear; close all;

data_folder = 'C:\Users\deu04\OneDrive\문서\MATLAB\Formation';
save_path = data_folder;
I_1C = 0.00382; %[A]

% MAT 파일 가져오기
slash = filesep;
files = dir([data_folder slash '*.mat']);

% formation 구조체 생성
num_files = length(files);
num_fields = 4;

formation(num_files) = struct('cycle', [], 'Q_chg', [], 'Q_dis', [], 'eff', []);

% 각 파일에 대해 q_chg 계산 및 formation 구조체에 저장
for file_index = 1:num_files
    % 선택한 파일 load
    fullpath_now = [data_folder slash files(file_index).name];
    load(fullpath_now);

    % Q값 계산하기
    for i = 1:length(data)
        data(i).Q = abs(trapz(data(i).t, data(i).I)) / 3600;
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

    formation(file_index).cycle = unique([data.cycle]);

    % formation 구조체 안에 Q_chg 지정
    for i = 1:length(step_chg)
        formation(file_index).Q_chg(i) = data(step_chg(i)).Q;
    end

    % formation 구조체 안에 Q_dis 지정
    for i = 1:length(step_dis)
        formation(file_index).Q_dis(i) = data(step_dis(i)).Q;
    end

    % formation 구조체 안에 eff 지정
    formation(file_index).eff = formation(file_index).Q_dis ./ formation(file_index).Q_chg;
end


% 평균값 구하기
num_cycles = length(formation(1).cycle); 

QC_mean = zeros(num_cycles, 1); 
QD_mean = zeros(num_cycles, 1); 

for i = 1:num_cycles
    QC_sum = 0; 
    QD_sum = 0;
    
    for j = 1:length(files)
        QC_sum = QC_sum + formation(j).Q_chg(i);
        QD_sum = QD_sum + formation(j).Q_dis(i);
    end
    
    QC_mean(i) = QC_sum / length(files);
    QD_mean(i) = QD_sum / length(files);
end

eff_mean = QD_mean ./ QC_mean;


% Error bar을 위한 min,max 구하기

min_cycle_chg = zeros(4, 1);
max_cycle_chg = zeros(4, 1);

for i = 1:4
    min_val_chg = inf;  % Set initial minimum value to infinity
    max_val_chg = -inf; % Set initial maximum value to negative infinity
    for j = 1:3
        min_val_chg = min(min_val_chg, formation(j).Q_chg(i));  % Update minimum value
        max_val_chg = max(max_val_chg, formation(j).Q_chg(i));  % Update maximum value
    end
    min_cycle_chg(i) = min_val_chg;  % Store the minimum value
    max_cycle_chg(i) = max_val_chg;  % Store the maximum value
end

min_cycle_dis = zeros(4, 1);
max_cycle_dis = zeros(4, 1);

for i = 1:4
    min_val_dis = inf;  % Set initial minimum value to infinity
    max_val_dis = -inf; % Set initial maximum value to negative infinity
    for j = 1:3
        min_val_dis = min(min_val_dis, formation(j).Q_dis(i));  % Update minimum value
        max_val_dis = max(max_val_dis, formation(j).Q_dis(i));  % Update maximum value
    end
    min_cycle_dis(i) = min_val_dis;  % Store the minimum value
    max_cycle_dis(i) = max_val_dis;  % Store the maximum value
end



% Plot
figure;

yyaxis left;
errorbar(1:num_cycles, QC_mean, QC_mean-min_cycle_chg, max_cycle_chg-QC_mean, 'b-o');
axis([0 4 0 0.006])

hold on;
plot(1:num_cycles, QD_mean, 'r-o');
errorbar(1:num_cycles, QD_mean, QD_mean-min_cycle_dis, max_cycle_dis-QD_mean, 'r');
axis([0 4 0 0.006])

ylabel('QC_{mean}, QD_{mean}');
yyaxis right;
plot(1:num_cycles, eff_mean, 'g-o');
axis([0 4 0 1.5])

ylabel('eff_{mean}');
xlabel('cycle\_index');
legend('QC', 'QD', 'eff', 'Location', 'Best'); 
