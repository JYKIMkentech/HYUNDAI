%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 코드 : MCT-1~6 데이터를 불러와 처리
%        1) OCV 시트 로드 및 중복 제거
%        2) MCT-1~6 각 시트에서 초기 Rest(휴지) 지점 탐색
%        3) MCT-1에 대해서 Coulomb Counting으로 SOC 그래프 추가
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
ocvCellVoltage = dataOCV.CellVoltage;
socOCV         = dataOCV.SOC_OCV;
[uCellVoltage, idxUnique] = unique(ocvCellVoltage);
uSocOCV = socOCV(idxUnique);

disp('=== OCV 시트 데이터 로드 및 중복 제거 완료 ===');

%% 3) 미리 셀/테이블 초기화
mctCellData = cell(6,1);   % 각 MCT 데이터 테이블을 저장할 셀 (6행 × 1열)

% OCVMCT 테이블: 6행 2열 (각 열은 double 타입, 열 이름은 "OCV_SoC", "BMS_SoC")
OCVMCT = table('Size',[6 2], ...
               'VariableTypes',{'double','double'}, ...
               'VariableNames',{'OCV_SoC','BMS_SoC'});

%% 4) MCT-1 ~ MCT-6 시트를 순회하며 처리
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
    
    % (5) 필요한 열 추출
    time         = dataMCT.Time_s;
    speed        = dataMCT.Velocity_kmh;
    batteryCurr  = dataMCT.Current_A;
    packVoltage  = dataMCT.PackVoltage_V;
    cellVoltMax  = dataMCT.CellVoltMax_V;
    
    socDecimal   = dataMCT.SOC_decimal;
    socInteger   = dataMCT.SOC_integer;
    socMCT       = socDecimal + socInteger;  % BMS에서 표시하는 SOC(%)
    
    % (6) "초기 휴지 상태" 마지막 인덱스 찾기
    idxRest = find(batteryCurr ~= 0, 1) - 1;
    if isempty(idxRest) || idxRest < 1
        disp(['[MCT-' num2str(mctNumber) '] 초기 휴지 구간이 없거나 예외상황.']);
        idxRest = NaN;
    end
    
    %% (7) 기본 그래프(2×2 서브플롯) 그리기
    figure('Name',['MCT-' num2str(mctNumber) ' Data Plots'],'NumberTitle','off');
    
    % (a) Time vs Speed
    subplot(2,2,1);
    plot(time, speed, 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Speed (km/h)');
    title('Time vs Speed');
    grid on;
    
    % (b) Time vs Current
    subplot(2,2,2);
    plot(time, batteryCurr, 'LineWidth', 1.2, 'Color','r');
    hold on;
    % 마지막 휴지 구간 표시 (idxRest가 유효하면)
    if ~isnan(idxRest)
        plot(time(idxRest), batteryCurr(idxRest), 'ro','LineWidth',2,'MarkerSize',6);
    end
    xlabel('Time (s)');
    ylabel('Current (A)');
    title('Time vs Current');
    grid on;
    
    % (c) Time vs Pack Voltage
    subplot(2,2,3);
    plot(time, packVoltage, 'LineWidth', 1.2, 'Color','b');
    xlabel('Time (s)');
    ylabel('Pack Voltage (V)');
    title('Time vs Pack Voltage');
    grid on;
    
    % (d) Time vs SOC
    subplot(2,2,4);
    plot(time, socMCT, 'LineWidth', 1.2, 'Color','k');
    xlabel('Time (s)');
    ylabel('SOC (%)');
    title('Time vs SOC');
    grid on;
    
    sgtitle(['MCT-' num2str(mctNumber) ' Driving Cycle Data'],...
            'FontWeight','bold','FontSize',12);
    
    %% (8) Rest 지점의 CellVoltMax_V를 이용한 SOC 역추적(OCV) 및 비교
    if ~isnan(idxRest)
        restVoltage = cellVoltMax(idxRest);
        
        % 보간으로 SOC 추정 (OCV 기반)
        socEstimate = interp1(uCellVoltage, uSocOCV, restVoltage, 'linear', 'extrap');
        socFromCSV  = socMCT(idxRest);
        
        % (i) 결과 콘솔 출력
        fprintf('\n[MCT-%d]\n', mctNumber);
        fprintf('  - 초기 휴지 구간 마지막 인덱스: %d\n', idxRest);
        fprintf('  - 해당 시점 CellVoltMax_V: %.4f V\n', restVoltage);
        fprintf('  - OCV 보간 SOC: %.2f %%\n', socEstimate);
        fprintf('  - CSV 기록 SOC: %.2f %%\n', socFromCSV);
        fprintf('  => 차이: %.2f %%p\n', socEstimate - socFromCSV);
        
        % (ii) OCVMCT 테이블에 저장
        OCVMCT.OCV_SoC(mctNumber) = socEstimate;  % 첫 번째 열: OCV 기반 SOC
        OCVMCT.BMS_SoC(mctNumber) = socFromCSV;   % 두 번째 열: BMS SOC
        
        % (iii) MCT-1의 경우 쿨롱 카운팅(Coulomb Counting) 계산 및 비교
        if mctNumber == 1
            % (a) 셀 정격 용량 (예시 값, 실제로 수정 필요)
            capacityAh = 56.2396;  
            
            % (b) 쿨롱카운팅용 SOC 벡터 초기화
            coulombSoc = zeros(size(time));
            
            % 초기 구간(1 ~ idxRest)까지 모두 OCV 기반 SOC로 설정
            coulombSoc(1 : idxRest) = socEstimate;
            
            % (c) 시간 간격 (초 단위)
            dt_sec = diff(time);
            
            % (d) idxRest+1부터 쿨롱 카운팅 적분
            for k = (idxRest+1) : length(time)
                dt = dt_sec(k-1);      % [s]
                I  = batteryCurr(k-1); % [A]
                
                % 방전 전류(+) → SOC 감소 (가정)
                dAh = I * (dt / 3600);          % [Ah]
                coulombSoc(k) = coulombSoc(k-1) - (dAh / capacityAh)*100;
            end
            
            % (e) 쿨롱 카운팅 결과 그래프
            figure('Name','MCT-1 SOC Comparison','NumberTitle','off');
            plot(time, socMCT,'LineWidth',1.2,'DisplayName','BMS SOC'); 
            hold on;
            plot(time, coulombSoc,'LineWidth',1.2,'DisplayName','Coulomb Counting SOC');
            xlabel('Time (s)');
            ylabel('SOC (%)');
            grid on; 
            legend('Location','best');
            title('SOC Comparison for MCT-1');
        end
        
    else
        % 휴지 구간이 없으면 NaN으로 저장
        OCVMCT.OCV_SoC(mctNumber) = NaN;
        OCVMCT.BMS_SoC(mctNumber) = NaN;
    end
    
end

disp('=== 모든 MCT 시트 데이터 처리 및 그래프 그리기 완료 ===');

%% (9) OCV 시트 그래프 (가로:세로 비율 2:1)
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

%% (10) 최종적으로 mctCellData, OCVMCT, OCV 관련 변수를 저장
save('MCT_Results.mat', ...
     'mctCellData', ...     % 1~6번 MCT의 전체 테이블들을 담은 cell
     'OCVMCT', ...          % Rest 시점 SOC 비교 결과 테이블
     'dataOCV', ...         % OCV 전체 원본 테이블 (필요하면)
     'socOCV', ...          % OCV용 SOC 벡터
     'ocvCellVoltage', ...  % OCV용 전압 벡터
     'uSocOCV', ...         % 중복 제거된 SOC
     'uCellVoltage');       % 중복 제거된 전압

disp('=== MCT_Results.mat 파일로 저장 완료 ===');


