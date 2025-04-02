close all
clear all
clc

% compare NCP btl and pigi september
% Ben Lowin
% October 25th 2024

%%  Load data

load PL02_NCP_01.mat


%% translate the ncp to /m^3
NCP_m3=pigi_dat.NCP./3.5;

%% plot figure

figure
plot(pigi_dat.datetime,NCP_m3)

ylabel(['Net Comunity Production' newline ' (mmol O_{2} / m^{3} / day)'])








