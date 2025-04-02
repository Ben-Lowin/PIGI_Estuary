

close all
clear all
clc

%make a new TDGP dataset
% Ben Lowin
% October 28th 2024

%%
load PL01_NCP_03.mat
PL01 = pigi_dat;

load PL02_NCP_03.mat
PL02 = pigi_dat;

load PL03_NCP_03.mat
PL03 = pigi_dat;

clear pigi_dat

% 300:52650
% 52651 :74562
% 74563:end
% point ranges for reference

%% plot the O2 to TDG plot, showing the strong linear relationship

O2_raw = PL02.raw.o2uM(300:52650);
Tdg_raw = PL02.raw.tp(300:52650);

x = [ones(size(O2_raw)), O2_raw];
y= Tdg_raw;

%regress needs x to have a column of ones

[b,bint,r,rint,stats] = regress(y,x);

mdl = fitlm(O2_raw,Tdg_raw);

figure
plot(mdl)
xlabel('Oxygen (umol/L)')
ylabel('Total disolved Gas (mbar)')
text(125,1010, ['R^2 = ' num2str(stats(1)) newline ...
    'y = 0.752 x + 846.94'])



%% atempt with O2 data as magnitude
%does not work, moving on
% O2_magnitude = O2_raw./mean(O2_raw,'omitnan');
% m_tdg = mean(Tdg_raw,'omitnan');
% 
% m_tdg = m_tdg.*O2_magnitude;


%%  Check the regression against raw O2 data

r_tdg = O2_raw.*b(2)+b(1);

figure
plot(Tdg_raw)
hold on
plot(r_tdg)
legend('raw', 'regressin')
hold off

%% giving the regression a try with unknown data

rf_tdg = pigi_dat.raw.o2uM.*b(2)+b(1);


figure
plot(pigi_dat.datetime,pigi_dat.raw.tp)
hold on
plot(pigi_dat.datetime,rf_tdg)
legend('raw', 'regressin')
hold off





