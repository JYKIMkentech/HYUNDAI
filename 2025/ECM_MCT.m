%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script Example
%   - 사용자에게 1~6 중 하나를 입력받아 해당 "MCT-#" 시트 데이터를 불러옴
%   - OCV 시트도 항상 불러와서 SOC vs (Cell 전압, Pack 전압) 그래프를 그림
%   - MCT 데이터는 기존과 같이 5행 헤더, 6행부터 실제 데이터라고 가정
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 0) 사용자 입력 받기
mctNumber = input('분석할 MCT 번호를 입력하세요 (1~6): ');  % 예: 1 -> MCT-1 시트
sheetNameMCT = ['MCT-' num2str(mctNumber)];

%% 1) 엑셀 파일 경로/이름 설정
filename = 'G:\공유 드라이브\BSL-Data\Data\Hyundai_dataset\현대차파우치셀 (rOCV,Crate)\NE_MCT25oC_HPPC25oC_OCV_KENTECH_송부.xlsx';

%% 2) MCT 시트 데이터 불러오기
optsMCT = detectImportOptions(filename, 'Sheet', sheetNameMCT, 'VariableNamingRule','preserve');
optsMCT.VariableNamesRange = 'A5:J5';  % 헤더(컬럼명)가 5행에 있다고 가정
optsMCT.DataRange          = 'A6';     % 데이터는 6행부터 시작한다고 가정
dataMCT                    = readtable(filename, optsMCT);

% 불러온 데이터의 열 이름 등을 확인(필요시)
disp('=== MCT 데이터 상위 5행 ===');
disp(head(dataMCT, 5));
disp('=== MCT 열(변수) 이름 목록 ===');
disp(dataMCT.Properties.VariableNames);

%% 3) MCT 시트에서 실제 사용할 열 이름 지정 및 데이터 추출
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

time        = dataMCT.Time_s;
speed       = dataMCT.Velocity_kmh;
batteryCurr = dataMCT.Current_A;
packVoltage = dataMCT.PackVoltage_V;

% SOC = (소수부 + 정수부)
socDecimal = dataMCT.SOC_decimal;
socInteger = dataMCT.SOC_integer;
socMCT     = socDecimal + socInteger;

%% 4) OCV 시트 데이터 불러오기 (항상 로드)
sheetNameOCV = 'OCV';  % 실제 파일의 OCV 시트 이름에 맞게 수정
optsOCV = detectImportOptions(filename, 'Sheet', sheetNameOCV, 'VariableNamingRule','preserve');
% 예: OCV 시트에서 1행에 헤더가 있고, 2행부터 데이터라고 가정
optsOCV.DataRange = 'A3';
dataOCV = readtable(filename, optsOCV);

% 사용 중인 예시 스크린샷에 따르면 열이 (SOC, Cell, BSA) 순으로 있다고 가정
% 실제 파일 형식에 맞게 열 이름을 수정
dataOCV.Properties.VariableNames{1} = 'SOC_OCV';     % SOC (%)
dataOCV.Properties.VariableNames{2} = 'CellVoltage'; % 단일 Cell 전압 (예: V)
dataOCV.Properties.VariableNames{3} = 'PackVoltage'; % Pack 전압 or BSA (예: V)

socOCV        = dataOCV.SOC_OCV;
ocvCellVoltage = dataOCV.CellVoltage;
ocvPackVoltage = dataOCV.PackVoltage;

%% 5) 그래프 그리기

% (1) MCT 데이터 서브플롯
figure('Name',['MCT-' num2str(mctNumber) ' Data Plots'],'NumberTitle','off');

subplot(2,2,1);
plot(time, speed, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Speed (km/h)');
title('Time vs Speed');
grid on;

subplot(2,2,2);
plot(time, batteryCurr, 'LineWidth', 1.2, 'Color','r');
xlabel('Time (s)');
ylabel('Current (A)');
title('Time vs Current');
grid on;

subplot(2,2,3);
plot(time, packVoltage, 'LineWidth', 1.2, 'Color','b');
xlabel('Time (s)');
ylabel('Pack Voltage (V)');
title('Time vs Pack Voltage');
grid on;

subplot(2,2,4);
plot(time, socMCT, 'LineWidth', 1.2, 'Color','k');
xlabel('Time (s)');
ylabel('SOC (%)');
title('Time vs SOC');
grid on;

sgtitle(['MCT-' num2str(mctNumber) ' Driving Cycle Data'],...
        'FontWeight','bold','FontSize',12);

% (2) OCV 데이터 그래프 (SOC vs CellVoltage, PackVoltage)
figure('Name','OCV Data','NumberTitle','off');

% Cell Voltage vs SOC
subplot(1,2,1);
plot(socOCV, ocvCellVoltage, 'o-','LineWidth',1.2);
xlabel('SOC (%)');
ylabel('Cell Voltage (V)');
title('SOC vs Cell OCV');
grid on;

% Pack Voltage vs SOC
subplot(1,2,2);
plot(socOCV, ocvPackVoltage, 's-','LineWidth',1.2,'Color','m');
xlabel('SOC (%)');
ylabel('Pack Voltage (V)');
title('SOC vs Pack OCV');
grid on;

sgtitle('OCV Test Data','FontWeight','bold','FontSize',12);

disp('--- 그래프가 모두 그려졌습니다. ---');

% 1) OCV 전압 벡터 중복 제거
[uCellVoltage, idx] = unique(ocvCellVoltage);

% 2) 첫 번째 행의 CellVoltMax_V를 기준으로 SOC를 보간
socEstimate = interp1(uCellVoltage, socOCV(idx), dataMCT{1,'CellVoltMax_V'}, 'linear','extrap');



