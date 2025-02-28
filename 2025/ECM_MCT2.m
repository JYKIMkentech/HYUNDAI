%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 코드 : MCT-1~6 데이터를 불러와 처리
%        1) OCV 시트 로드 및 중복 제거
%        2) MCT-1~6 각 시트에서 초기 Rest(휴지) 지점 탐색
%        3) Coulomb Counting(C.C.) SOC 계산 및 BMS SOC와 비교
%        4) 결과 mctCellData, OCVMCT, etc. 저장
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1) 엑셀 파일 경로/이름 설정
filename = 'G:\공유 드라이브\BSL-Data\Data\Hyundai_dataset\현대차파우치셀 (rOCV,Crate)\NE_MCT25oC_HPPC25oC_OCV_KENTECH_송부.xlsx';

%% 2) OCV 시트 데이터 불러오기 (1회만)
sheetNameOCV = 'OCV';  % 실제 파일의 OCV 시트 이름
optsOCV = detectImportOptions(filename, 'Sheet', sheetNameOCV, 'VariableNamingRule','preserve');
optsOCV.DataRange = 'A3';
dataOCV = readtable(filename, optsOCV);

% 열 이름 지정 (실제 파일 열 순서/내용 확인 후 수정)
dataOCV.Properties.VariableNames{1} = 'SOC_OCV';
dataOCV.Properties.VariableNames{2} = 'CellVoltage';
dataOCV.Properties.VariableNames{3} = 'PackVoltage';

% OCV 중복 제거
ocvCellVoltage = dataOCV.CellVoltage;   % 셀 기준 전압
socOCV         = dataOCV.SOC_OCV;
[uCellVoltage, idxUnique] = unique(ocvCellVoltage);
uSocOCV = socOCV(idxUnique);

disp('=== OCV 시트 데이터 로드 및 중복 제거 완료 ===');


%% 3) 미리 변수/테이블 초기화
mctCellData = cell(6,1);   % 각 MCT 데이터 테이블을 저장할 셀 (6행 × 1열)

% OCVMCT 테이블(초기 Rest 시점 SOC 비교용): 6행 2열
OCVMCT = table('Size',[6 2], ...
               'VariableTypes',{'double','double'}, ...
               'VariableNames',{'OCV_SoC','BMS_SoC'});

           
%% (추가) 배터리 정보(예시)
%  - 직렬: 192개, 병렬: 2개
numSeries   = 192;
numParallel = 2;
Q_batt      = 56.2396;  % [Ah]


%% 4) MCT-1 ~ MCT-6 시트 순회하며 처리
for mctNumber = 1:6
    
    % (1) 시트명 설정
    sheetNameMCT = ['MCT-' num2str(mctNumber)];
    
    % (2) 시트 데이터 불러오기
    optsMCT = detectImportOptions(filename, 'Sheet', sheetNameMCT, 'VariableNamingRule','preserve');
    optsMCT.VariableNamesRange = 'A5:J5';  % 예: 5행에 헤더(컬럼명)
    optsMCT.DataRange          = 'A6';     % 예: 6행부터 실제 데이터
    dataMCT = readtable(filename, optsMCT);
    
    % (3) 열 이름 설정 (실제 파일 포맷에 맞게 수정)
    dataMCT.Properties.VariableNames{1}  = 'Time_s';
    dataMCT.Properties.VariableNames{2}  = 'Velocity_kmh';
    dataMCT.Properties.VariableNames{3}  = 'Current_A';
    dataMCT.Properties.VariableNames{4}  = 'PackVoltage_V';
    dataMCT.Properties.VariableNames{5}  = 'CellVoltMax_V';
    dataMCT.Properties.VariableNames{6}  = 'TempMax';
    dataMCT.Properties.VariableNames{7}  = 'CellVoltMin_V';
    dataMCT.Properties.VariableNames{8}  = 'TempMin';
    dataMCT.Properties.VariableNames{9}  = 'SOC_decimal';
    dataMCT.Properties.VariableNames{10} = 'SOC_integer';
    
    % (4) mctCellData에 저장
    mctCellData{mctNumber} = dataMCT;  % 각 MCT 테이블을 셀로 저장
    
    %% (5) 필요한 열 추출
    time_s      = dataMCT.Time_s;
    speed_kmh   = dataMCT.Velocity_kmh;
    packCurrent = dataMCT.Current_A;     % (양수: 방전)
    packVoltage = dataMCT.PackVoltage_V;
    
    socDecimal  = dataMCT.SOC_decimal;
    socInteger  = dataMCT.SOC_integer;
    SOC_bms     = socDecimal + socInteger;  % [%] BMS SOC
    
    %% (6) 초기 휴지 상태(전류=0) 마지막 인덱스 찾기
    idxRest = find(packCurrent ~= 0, 1) - 1;
    if isempty(idxRest) || idxRest < 1
        disp(['[MCT-' num2str(mctNumber) '] 초기 휴지 구간이 없거나 예외상황.']);
        idxRest = NaN;
    end
    
    %% (7) Coulomb Counting용 변수 계산
    %  - Pack → Cell 전압(192 직렬), Pack → Cell 전류(2 병렬)
    cellVoltage_meas = packVoltage / numSeries;
    cellCurrent       = packCurrent / numParallel;   % [A]
    
    % 시간 간격
    dt = [0; diff(time_s)];
    
    % 7.1) 초기 SOC 계산 (OCV 기반)
    %      "처음으로 전류가 비영이 되는 시점" 바로 이전 전압을 초기 전압으로 사용
    idx_firstNonZero = find(cellCurrent ~= 0, 1, 'first');
    if isempty(idx_firstNonZero)
        idx_init = 1;
        warning('전류가 전체 구간에서 0입니다. 초기 인덱스를 1로 설정합니다.');
    else
        idx_init = idx_firstNonZero - 1;
        if idx_init < 1
            idx_init = 1;
            warning('첫 지점부터 전류가 0이 아니므로, 초기 인덱스를 1로 설정합니다.');
        end
    end
    
    cellVoltage_init = cellVoltage_meas(idx_init);
    % 중복 제거된 OCV 데이터를 이용하여 SOC0 보간
    SOC0 = interp1(uCellVoltage, uSocOCV, cellVoltage_init, 'linear','extrap');
    
    % 7.2) 누적 전류 적분 → SOC(t)
    %      cumtrapz(time, current)를 이용하면 적분값(A·s)이 나오므로,
    %      Ah로 변환 후 정격 용량 대비 백분율을 구한다.
    charge_integral = cumtrapz(time_s, cellCurrent);  % [A·s]
    SOC_cc = SOC0 - (charge_integral/(Q_batt*3600))*100;  % [%]
    
    
    %% (8) 2×2 서브플롯으로 데이터 시각화
    figure('Name',['MCT-' num2str(mctNumber) ' Data Plots'],'NumberTitle','off');
    
    % (a) Time vs Speed
    subplot(2,2,1);
    plot(time_s, speed_kmh, 'LineWidth', 1.2, 'Color','b');
    xlabel('Time (s)');
    ylabel('Speed (km/h)');
    title('Time vs Speed');
    grid on;
    
    % (b) Time vs Current
    subplot(2,2,2);
    plot(time_s, packCurrent, 'LineWidth', 1.2, 'Color','r');
    hold on;
    % 초기 휴지 구간 마지막 인덱스 표시
    if ~isnan(idxRest)
        plot(time_s(idxRest), packCurrent(idxRest), 'ro','LineWidth',2,'MarkerSize',6);
    end
    xlabel('Time (s)');
    ylabel('Current (A)');
    title('Time vs Current');
    grid on;
    
    % (c) Time vs Pack Voltage
    subplot(2,2,3);
    plot(time_s, packVoltage, 'LineWidth', 1.2, 'Color','b');
    xlabel('Time (s)');
    ylabel('Pack Voltage (V)');
    title('Time vs Pack Voltage');
    grid on;
    
    % (d) Time vs SOC (BMS SOC vs CC SOC)
    subplot(2,2,4);
    hold on;
    plot(time_s, SOC_bms, 'LineWidth', 1.3, 'Color','b', 'DisplayName','BMS SOC');
    plot(time_s, SOC_cc,  'LineWidth', 1.3, 'Color','r', 'DisplayName','CC SOC');
    xlabel('Time (s)');
    ylabel('SOC (%)');
    title('Time vs SOC');
    legend('Location','best');
    box on;
    grid on;
    
    sgtitle(['MCT-' num2str(mctNumber) ' Driving Cycle Data'],...
            'FontWeight','bold','FontSize',12);
    
    %% (9) 초기 휴지 상태의 CellVoltMax_V를 이용한 SOC 역추적(OCV) → OCVMCT 기록
    if ~isnan(idxRest)
        restVoltage = dataMCT.CellVoltMax_V(idxRest);
        
        % 셀 기준 OCV 중복 제거 데이터를 이용
        socEstimate = interp1(uCellVoltage, uSocOCV, restVoltage, 'linear', 'extrap');
        socFromCSV  = SOC_bms(idxRest);  % BMS 기록 SOC
        
        % 콘솔 출력
        fprintf('\n[MCT-%d]\n', mctNumber);
        fprintf('  - 초기 휴지 구간 마지막 인덱스: %d\n', idxRest);
        fprintf('  - 해당 시점 CellVoltMax_V: %.4f V\n', restVoltage);
        fprintf('  - OCV 보간 SOC: %.2f %%\n', socEstimate);
        fprintf('  - BMS SOC: %.2f %%\n', socFromCSV);
        fprintf('  => 차이: %.2f %%p\n', socEstimate - socFromCSV);
        
        % OCVMCT 테이블에 저장
        OCVMCT.OCV_SoC(mctNumber) = socEstimate;
        OCVMCT.BMS_SoC(mctNumber) = socFromCSV;
    else
        % 휴지 구간이 없으면 NaN
        OCVMCT.OCV_SoC(mctNumber) = NaN;
        OCVMCT.BMS_SoC(mctNumber) = NaN;
    end
    
end

disp('=== 모든 MCT 시트 데이터 처리 및 그래프 그리기 완료 ===');


%% (10) OCV 시트 그래프 (가로:세로 비율 2:1)
figure('Name','OCV Data','NumberTitle','off',...
       'Position',[100,100,1200,600]);  % (왼쪽위 X, 왼쪽위 Y, 폭, 높이)

% (a) Cell Voltage vs SOC
subplot(1,2,1);
plot(socOCV, ocvCellVoltage, 'o-','LineWidth',1.2);
xlabel('SOC (%)');
ylabel('Cell Voltage (V)');
title('SOC vs Cell OCV');
grid on;

% (b) Pack Voltage vs SOC
subplot(1,2,2);
plot(socOCV, dataOCV.PackVoltage, 's-','LineWidth',1.2,'Color','m');
xlabel('SOC (%)');
ylabel('Pack Voltage (V)');
title('SOC vs Pack OCV');
grid on;

sgtitle('OCV Test Data','FontWeight','bold','FontSize',12);

disp('--- 모든 그래프가 표시되었습니다. ---');


%% (11) 최종적으로 mctCellData, OCVMCT, OCV 관련 변수를 저장
save('MCT_Results.mat', ...
     'mctCellData', ...      % 1~6번 MCT의 전체 테이블들을 담은 cell
     'OCVMCT', ...           % Rest 시점 SOC 비교 결과 테이블
     'dataOCV', ...          % OCV 전체 원본 테이블 (필요하면)
     'socOCV', ...           % OCV용 SOC 벡터
     'ocvCellVoltage', ...   % OCV용 전압 벡터
     'uSocOCV', ...          % 중복 제거된 SOC
     'uCellVoltage');        % 중복 제거된 전압

disp('=== MCT_Results.mat 파일로 저장 완료 ===');


