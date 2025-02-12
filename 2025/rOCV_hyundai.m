clc;clear;close all

load('G:\공유 드라이브\BSL-Data\Data\Hyundai_dataset\현대차파우치셀 (rOCV,Crate)\rOCV(C 20)\NE_SOC-OCV_0.05C_Raw data.mat')

all_t = vertcat(data(12:84).t);
all_I = vertcat(data(12:84).I);

capacity_Ah = trapz(all_t, all_I) / 3600;

fprintf('총 용량 = %.4f Ah\n', capacity_Ah);

figure(1);
plot(vertcat(data.t), vertcat(data.I))
ylabel('Current')
hold on
yyaxis right
plot(vertcat(data.t), vertcat(data.V))
ylabel('Voltage')


figure(2);
plot(vertcat(data(11:84).t), vertcat(data(11:84).I))
ylabel('Current')
hold on
yyaxis right
plot(vertcat(data(11:84).t), vertcat(data(11:84).V))
ylabel('Voltage')

