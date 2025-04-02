%% ABOUT

%------------------------------------------------------------------------
% SCRIPT: process_pigi
% 
% ABOUT: 
%   Main script to compile and process PIGI data output from LabVIEW
%   acquisitions.
% 
%   Run the sections below to produce saved data structures and output
%   figures. NOTE: data should also be QC'd and calibrated before
%   performing calculations (section 4).
% 
% R. Izett
% rizett@eoas.ubc.ca
% Last modified: June 2020
%------------------------------------------------------------------------
        
%% SECTION 1: EXTRACT AND COMPILE ALL DATA FROM A SINGLE DEPLOYMENT

clearvars; close all; clc

%-------------------------------------------------------------------------
% INPUT INFORMATION

% Directory where PIGI processing code is saved
    code_dir = 'C:\Users\BenL\Documents\MATLAB\PIGI\Scripts';

% Directory where data are saved
    data_dir = 'C:\Users\BenL\Documents\MATLAB\PIGI\Dock_NCP_03\PIGI_deployment_data';
        
% Directory where output/processed data will be saved
    save_dir = 'C:\Users\BenL\Documents\MATLAB\PIGI\Dock_NCP_03\PIGI_out';
    
% Directory where output/processing data will be saved
    fig_dir = 'C:\Users\BenL\Documents\MATLAB\PIGI\Dock_NCP_03\PIGI_out';

% (Useful) Deployment information
    %REQUIRED:
    dep_info.name = 'PL_03_deployment'; %NOTE - do not use dashes in the name
    dep_info.year = 2024;
    dep_info.start = datestr(datenum(2024,12,08)); %datestring format (e.g. = datestr(datenum(2020,06,01));)
    dep_info.finish = datestr(datenum(2024,12,14)); %datestring format
    %OPTIONAL
    dep_info.id = 'PL03';
    dep_info.location = 'Preist Landing - Floating dock';
    dep_info.owner = 'Lowin';
    dep_info.collector = 'Lowin';
%-------------------------------------------------------------------------

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%
disp('--------------------------------------------')
disp('*START OF SECTION 1*')

%--- CD to directory where data are saved and extract raw data
    addpath(genpath(code_dir)); %add code directory to current path
    cd(data_dir)
    
    %--- Run function to extract data from PIGI file
        [pigi_dat,comnt,comnt_t] = pigi_processing(cd,dep_info); 

%     %--- Change the O2 % to %/100
%         pigi_dat.raw.o2sat=pigi_dat.raw.o2sat./100;
%         
    %--- Save data comments
        pigi_dat.comments.acquisition_comments = comnt;
        pigi_dat.comments.acquisition_comments_time = comnt_t;
        pigi_dat.comments.v1 = ['acquisition data extracted & compiled (' datestr(datenum(date),'yyyy-mmm-dd'),')'];
        clear comnt comnt_t
        disp('PIGI data extracted as: ')
        pigi_dat
        
    %--- Save RAW data
        sname = [dep_info.name,'_PIGI_RAW'];
        [sname] = save_version(save_dir,sname,pigi_dat,'pigi_dat');
        disp(['RAW data saved to: ' sname])
        clear sname


     %--- QA/QC for pigi data tp
        pigi_dat.raw.tp_fix = 0.752.* pigi_dat.raw.o2uM + 846.94;

% --- Plot original/raw data
    close all
    set(groot, 'defaultFigurePosition',[100 100 560 600]) %set figure position
    f=figure(1);

    %O2 saturation and GTD total pressure normalized to 1-atm (1013.25 mbar)
        subplot(3,1,1); hold on
        [ax,p1,p2] = plotyy(pigi_dat.time,pigi_dat.raw.o2sat,pigi_dat.time,100*pigi_dat.raw.tp./1013.25);
        set(p1,'color','k');set(p2,'color','r');
        set(ax(1),'ycolor','k');set(ax(2),'ycolor','r');
        format_plot(ax(1))
        ylabel(ax(1),'O2 percent-saturation [%]'); ylabel(ax(2),{'Normalized GTD total pressure'; '[mbar/1013.25]'});
        datetick(ax(1),'x');datetick(ax(2),'x');
        title('RAW PIGI DATA')

    %Optode temperature 
    %(a good indicator of when flow rate has stopped)
        ax(3)=subplot(3,1,2); hold on
        plot(pigi_dat.time,pigi_dat.opt_T,'k');
        format_plot
        ylabel('opt T [deg-C]');
        datetick('x');    

    %Flow rates
        ax(4) = subplot(3,1,3); hold on
        %note: flow rate data smoothed to reduce noise
        plot(pigi_dat.time,nanmoving_average(pigi_dat.prim_flow,20),'g')
        plot(pigi_dat.time,nanmoving_average(pigi_dat.inst_flow,20),'color',[.5 .5 .5]);
        format_plot
        ylabel('Flow rate [L/min]');
        datetick('x');
        legend('Primary chamber','Inst. loop','location','best')        

    linkaxes(ax,'x')    

    %--- Save figure
        sname = [dep_info.name,'_PIGI_RAW'];
        sname = [fig_dir,'/',sname,'_v',datestr(date,'yyyymmdd'),'.tif'];
        saveas(gcf,sname)
        disp(['RAW data figure saved to: ' sname])
        clear sname
    
disp('SECTION 1 COMPLETE: Data extracted and saved')
disp('--------------------------------------------')
disp(' ')

%% SECTION 2: OPTODE SALINITY COMPENSATION

clearvars -except data_dir save_dir fig_dir pigi_dat dep_info

%-------------------------------------------------------------------------
% LOAD TSG DATA AND OBTAIN TEMPERATURE, SALINITY AND TIME ARRAYS
    %NOTE: modify this section according to your T-S data source. If you do
        %not have T-S data available, set the values below to []. In this
        %case, the Optode data WILL NOT be T/S-corrected.
    %NOTE: Set ts_s to a single value (e.g. ts_s = 32) to perform
        %S-compensation at a single salinity value. Best-practice is to use
        %underway S, but this should provide an approximation in waters 
        %with approx. constant salinity. 
    %NOTE: you will need time, temperature and salinity arrays to perform
        %the calculations below.
    load("PL03_master.mat")
    ts_time = Inported_data.time';
    ts_t = Inported_data.temperature';
    ts_s = Inported_data.salinity';
    
% OPTIONAL: User-specified sea level pressure (SLP, mbar) data. Leave
% blank if you wish to use default (1020 mbar) value.
    slp_time = Inported_data.time';
    slp_mbar = Inported_data.barometeric';

% OPTIONAL : Load GPS data
%no gps, as it was stationary
      
% IDENTIFY FILE CONTAINING OPTODE SETTINGS 
    %NOTE: this file is an output file from the LabVIEW acquisition and
        %should be in the same location as the raw/underway data (i.e.
        %data_dir)
    cal.opt_settings = 'rv_primary_test_optode-settings.txt';
    %If the settins file is saved in the same location as the raw data, use
    %the following:
        %cal.opt_settings = [data_dir,'/optode-settings.txt'];
        
% IDENTIFY THE FUNCTION TO PERFORM OPTODE RE-CALCULATION
    func = 'svu'; 
    % 'svu' = Stern-Volmer Uchida (Uchida et al., 2008)
    %   (*NOTE: requires that Optode foil coefficients have been 
    %   obtained using a multi-point calibration)
    % 'ssv' = Simplified Stern-Volmer (GEOMAR/Bittig et al.,2018); 
    %   (*NOTE: requires that Optode fiol coefficients have been derivd 
    %   using this formula during mulit-point  calibration)
    % 'poly' = emperical polynomial 
    %   (Use this function if no foil coefficients are provided / if
    %   the sensor has not been multi-point calibrated)
%-------------------------------------------------------------------------

pigi_dat.datetime = datetime(pigi_dat.time, 'ConvertFrom', 'datenum');

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%
disp('--------------------------------------------')
disp('*START OF SECTION 2*')

%--- Interpolate SLP data to match PIGI data
    if ~isempty(slp_time)
        uw_slp = interp1(slp_time, slp_mbar, pigi_dat.datetime);
    else
        uw_slp = repmat(1020, size(pigi_dat.time));
        %update to change the assumed slp
    end
    
%--- Inpterpolate temperature & salinity data to match PIGI data
    %IF no T/S data provided, use Opt. T and 0 sal. values
    if isempty(ts_s)
        uw_tem = pigi_dat.opt_T;
        uw_sal = repmat(0,size(pigi_dat.time));
    %IF sal. estimated with single value, use that value and set SST to
    %Opt. T
    elseif numel(ts_s) == 1
        uw_tem = pigi_dat.opt_T;
        uw_sal = repmat(ts_s,size(pigi_dat.time));
    %IF underway T/S data provided, interpolated data to PIGI data
    %resolution
    else
        uw_tem = interp1(ts_time, ts_t, pigi_dat.datetime);
        uw_sal = interp1(ts_time, ts_s, pigi_dat.datetime);
    end
   
    %--- Perform Optode temperature & salinity compensation and re-calculate O2-conc. & O2-sat    
    % per Bittig et al. (2015; 2018), Uchida et al. (2008; 2010) & Optode
    % Manual
    
    % Use "best" temperature and salinity data to perform T and S
    % compensation: Temperature readings should be as close to the Optode / 
    % PIGI system as possible. TSG SST data are  preferred over Optode T
    % sensor because the TSG sensor will be more accurate and have a 
    % quicker response time.
    %IF no underway S data provided, copy raw data
    if all(uw_sal ==0)
        pigi_dat.ts_cor.o2uM = pigi_dat.raw.o2uM;
        pigi_dat.ts_cor.o2sat = pigi_dat.raw.o2sat;
        pigi_dat.ts_cor.po2 = O2ctoO2p(pigi_dat.ts_cor.o2uM,uw_tem,uw_sal,0);       
        
        disp('Optode data NOT T/S-corrected');
    %IF underway S provided (or estimated), peform S and T compensation
    else
        [pigi_dat.ts_cor.o2uM,pigi_dat.ts_cor.o2sat,pigi_dat.ts_cor.po2] = ...
            optode_recalc(pigi_dat.raw.opt_calPhase,uw_tem,uw_sal,uw_slp,cal.opt_settings,func);
        
        disp('Optode data T/S-corrected');        
    end
    pigi_dat.ts_cor.sal = uw_sal;
    pigi_dat.ts_cor.temp = uw_tem;
    pigi_dat.air.slp = uw_slp;
    pigi_dat.cal = cal;

    pigi_dat.ts_cor.windspeed = interp1(ts_time, Inported_data.wind_speed', pigi_dat.datetime);
    pigi_dat.ts_cor.air_temp = interp1(ts_time, Inported_data.air_temp', pigi_dat.datetime);
    pigi_dat.ts_cor.pH = interp1(ts_time, Inported_data.pH', pigi_dat.datetime);
    pigi_dat.ts_cor.density = interp1(ts_time, Inported_data.density', pigi_dat.datetime);
    pigi_dat.ts_cor.DO_hanna = interp1(ts_time, Inported_data.disolved_oxygen', pigi_dat.datetime);
    pigi_dat.ts_cor.sun_hours = interp1(ts_time, Inported_data.sun_hours', pigi_dat.datetime);
    pigi_dat.ts_cor.solar_rad = interp1(ts_time, Inported_data.solar_rad', pigi_dat.datetime);

    
%--- Save data 
    pigi_dat.comments.v2 = ['Optode data T/S-corrected (' datestr(datenum(date),'yyyy-mmm-dd'),')'];
    sname = [dep_info.name,'_PIGI_RAW'];
    [sname] = save_version(save_dir,sname,pigi_dat,'pigi_dat');
    disp(['RAW data (with Opt. S-compensation) saved to: ' sname])
    clear sname
    
%--- Plot raw and S-compensated O2
    figure;
    sp(1)=subplot(3,1,1); hold on; %O2-conc
        plot(pigi_dat.time,pigi_dat.ts_cor.o2uM,'c.-');
        plot(pigi_dat.time,pigi_dat.raw.o2uM,'k');
        format_plot
        ylabel('O2 concentration [\muM]');
        datetick('x');    
        legend('T/S-corrected','Raw','orientation','horizontal','location','best')
        title('T/S-Corrected PIGI DATA')
    
    sp(2)=subplot(3,1,2); hold on; %O2-saturation
        plot(pigi_dat.time,pigi_dat.ts_cor.o2sat,'c.-');
        plot(pigi_dat.time,pigi_dat.raw.o2sat,'k');        
        format_plot
        ylabel('O2 saturation state [%]');
        datetick('x');    
        
    sp(3)=subplot(3,1,3); hold on; %Salinity
        plot(pigi_dat.time,pigi_dat.ts_cor.sal,'k')
        format_plot
        ylabel('Salinity [PSU]')
        datetick('x')
    
    linkaxes(sp,'x')    

    %--- Save figure
        sname = [dep_info.name,'_PIGI_O2-TS-Comp'];
        sname = [fig_dir,'/',sname,'_v',datestr(date,'yyyymmdd'),'.tif'];
        saveas(gcf,sname)
        disp(['O2 T/S-compensation figure saved to: ' sname])
        clear sname

disp('SECTION 2 COMPLETE: Optode O2 S-compensated')
disp('--------------------------------------------')
disp(' ')
           
%% SECTION 3: DATA QC'ing and O2 calibration

clearvars -except data_dir save_dir fig_dir pigi_dat dep_info

%-------------------------------------------------------------------------
% Data should be QC'd (de-spiked, smoothed etc.) and O2
% and total dissolved gas pressure should be calibrated BEFORE 
% calculating N2 from GTD pressure. 

% This section has been left blank for the user to populate. 
% See Izett et al. 2020 for recommendations data processing steps.
%-------------------------------------------------------------------------

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%
disp('--------------------------------------------')
disp('*START OF SECTION 3*')

%--- QC data: filter, de-spike etc.
 %[x,y]=ginput;

%Fid where temp = 0

bad_data = find(pigi_dat.opt_T ==0);

%set bad data to NaN
pigi_dat.ts_cor.o2uM(bad_data)=nan;
pigi_dat.ts_cor.o2sat(bad_data)=nan;
pigi_dat.ts_cor.po2(bad_data)=nan;
pigi_dat.ts_cor.o2sat(bad_data)=nan;
pigi_dat.raw.tp(bad_data)=nan;
pigi_dat.raw.tp_fix(bad_data)=nan;
    
%--- Optode calibtation
%based on Winkler titrations - fix the offset
pigi_dat.qc_cal.o2uM    =1.1638 .* pigi_dat.ts_cor.o2uM - 2.0634;

pigi_dat.qc_cal.o2sat = O2ctoO2s(pigi_dat.qc_cal.o2uM , pigi_dat.ts_cor.temp, ...
    pigi_dat.ts_cor.sal,0,pigi_dat.air.slp);%saturation (%/100)

pigi_dat.qc_cal.po2 = O2ctoO2p(pigi_dat.qc_cal.o2uM,pigi_dat.ts_cor.temp, ...
    pigi_dat.ts_cor.sal,0); %partial pressure (mbar)

%--- Calibrate GTD total dissolved gas pressure data for offset using
    %atmospheric / in-air measurements.

%--- Correct GTD total dissolved gas pressure data 
    %for partial pressure increase associated with SST warming in seawater
    %intake lines.

%--- Save data 
    %FOR NOW: copy t/s-cor. data to qc/cal data
    %MODIFY below to save QC'd data. 

    pigi_dat.qc_cal.sal     = pigi_dat.ts_cor.sal;
    pigi_dat.qc_cal.temp    = pigi_dat.ts_cor.temp;
    pigi_dat.qc_cal.tp      = pigi_dat.raw.tp;
    pigi_dat.qc_cal.tp_fix      = pigi_dat.raw.tp_fix;
    
    pigi_dat.comments.v3 = ['Data cleaned, filterd and O2 calibrated (' datestr(datenum(date),'yyyy-mmm-dd'),')'];
    
    sname = [dep_info.name,'_PIGI_RAW'];
    [sname] = save_version(save_dir,sname,pigi_dat,'pigi_dat');
    disp(['QC''d ata saved to: ' sname])
    clear sname
    
disp('SECTION 3 COMPLETE: Data cleaned, filtered & O2 calibrated')
disp('--------------------------------------------')
disp(' ')    

%% SECTION 4: CALCULATE N2 AND delta-O2/N2 

clearvars -except data_dir save_dir fig_dir pigi_dat dep_info cal

%-------------------------------------------------------------------------
% SPECIFY DATA TO USE IN N2 CALCULATIONS. QC'd and calibrated data
% recommended (see section above)
    input_dat.o2  = pigi_dat.qc_cal.o2sat; %O2 %/100
    input_dat.tp  = pigi_dat.qc_cal.tp_fix; %total dissolved gas pressure (mabr)
    input_dat.sst = pigi_dat.qc_cal.temp; %SST (C)
    input_dat.sal = pigi_dat.qc_cal.sal; %Salinity (PSU)
%-------------------------------------------------------------------------

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%
disp('--------------------------------------------')
disp('*START OF SECTION 4*')
       
%--- Calculate N2 and dO2/N2
    %--- Input data
        %Confirm O2 in correct units
            if nanmean(input_dat.o2) > 10
                input_dat.o2 = input_dat.o2/100;
            end
        
        %Make sure all data are same size
            air.slp = pigi_dat.air.slp;
            
            a_fds = fields(air);
            d_fds = fields(input_dat);
            for kk = 1:numel(a_fds)
                air.(a_fds{kk}) = reshape(air.(a_fds{kk}),size(input_dat.(d_fds{1})));
            end
            for kk = 1:numel(d_fds)
                input_dat.(d_fds{kk}) = reshape(input_dat.(d_fds{kk}),size(input_dat.(d_fds{1})));
            end
            clear kk a_fds d_fds
        
        %Set Argon saturation to [] (unless Ar or O2/Ar data are available)
        ar = [];
        
        %Set pCO2 to [] (unless it is known)
        pco2 = [];

    %--- Perform calculations
        [dat,p,air] = n2_from_tp(input_dat.tp,input_dat.o2,ar,pco2,input_dat.sst,input_dat.sal,air);
        
        pigi_dat.n2_calc = dat;
        pigi_dat.n2_calc.o2sat = input_dat.o2;
        pigi_dat.n2_calc.tp = input_dat.tp;
        pigi_dat.n2_calc.sst = input_dat.sst;
        pigi_dat.n2_calc.sal = input_dat.sal;
        pigi_dat.n2_calc.gas_part_pres = p;
        
        pigi_dat.air = air;
        clear dat p
        disp('N2 calculated');
        
% --- Save data 
    %Re-order data fields
    pigi_dat = orderfields(pigi_dat, {'time','datetime','gps_time','gps_lat','gps_long','opt_T','gtd_T','inst_flow','prim_flow','flow_status',...
        'raw','ts_cor','qc_cal','n2_calc','air','cal','units','comments'})

    pigi_dat.comments.v4 = ['N2 and delO2/N2 calculated from calibrated O2 (' datestr(datenum(date),'yyyy-mmm-dd'),')'];
    sname = [dep_info.name,'_PIGI_RAW'];
    [sname] = save_version(save_dir,sname,pigi_dat,'pigi_dat');
    disp(['RAW data (with N2-calculation) saved to: ' sname])
    clear sname
    
%--- O2 and N2 plots
    figure
    ax(1)=subplot(3,1,1); hold on; %O2 and N2 saturation
        plot(pigi_dat.time,100*(pigi_dat.n2_calc.o2sat-1),'b');
        plot(pigi_dat.time,100*(pigi_dat.n2_calc.n2sat-1),'r');
        format_plot
        ylabel('O2 and N2 supersaturation [%]')
        datetick('x');  
        legend('O_2','N_2')
        title('N2 CALCULATIONS')
        
    ax(2)=subplot(3,1,2); %delta-O2/N2
        plot(pigi_dat.time,pigi_dat.n2_calc.do2n2,'k');
        format_plot
        datetick('x');  
        ylabel('\DeltaO_2/N_2 [%]')
        
    ax(3)=subplot(3,1,3);
        plot(pigi_dat.time,pigi_dat.air.slp,'k')
        format_plot
        datetick('x');  
        ylabel('SLP (mbar)')
        
    linkaxes(ax,'x')    

    %--- Save figure
        sname = [dep_info.name,'_PIGI_N2_calc.mat'];
        sname = [fig_dir,'/',sname,'_v',datestr(date,'yyyymmdd'),'.tif'];
        saveas(gcf,sname)
        disp(['N2-calculations figure saved to: ' sname])
        clear sname
                
disp('SECTION 4 COMPLETE: N2 and delta-O2/N2 calculated')
disp('--------------------------------------------')
disp(' ')    
    