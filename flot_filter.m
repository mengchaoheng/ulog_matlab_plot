clear all;
close all;
clc;
addpath(genpath(pwd));
% you can run on terminal 
% ulog2csv log_8_2021-5-20-11-52-08.ulg 
% to get csv files
% =====================1==========================
% Install pyulog using pip first.https://github.com/PX4/pyulog.
% in MacOS, it maybe have been installed by the px4-dev
% =====================2==========================
% Make sure it has installed ulog2csv correctly (check the output of which ulog2csv in Linux/MacOS or where ulog2csv in Windows).
% =====================3==========================
% Change the following line in ulogviewver.m:
% command = ['!/usr/local/bin/ulog2csv ' ulgFileName '.ulg'];
% to 
% command = ['!your ulog2csv path' ulgFileName '.ulg'];
% and 
% ulgFileName = '00_41_22'; 
% to 
% ulgFileName = 'your log name'; 

% ----fig size, you have to change it for your fig

d2r=pi/180;
r2d=180/pi;
%% new method: defualt pid for endurance load wp , remove I for wind and max vel
% endurance:(17_14_37 1127s) 12_07_46 1137s   11_10_22 1140s  11_37_27 1120s  
% load:06_28_53 06_45_38  06_54_35
% wp:17_09_30 
% wind: 11_39_14  15_57_58
 % mav vel: 09_19_16 10_27_37

%% 2023 flight data
% Users/mch/Documents/FlightLog/log/2023-02-24/06_12_56.ulg
% 06_22_46
% 06_44_14
% 07_04_47
% 07_26_57

% 03_43_50
% /Users/mch/Documents/FlightLog/log/2023-02-23/07_44_32.ulg
% 07_45_53
% 07_48_07
% 07_50_17
% 08_06_35
% 08_29_05
% 08_39_27

% 20240614 08_15_39 08_39_28 09_05_41 09_15_41 
% 20240616  05_31_54  14_03_42 06_20_57 15_12_39
 %% pca test:16_37_44 12_15_41 15_03_31  df-1.12.3 (15_10_14)

 % real flight 07_18_37  12_03_09   (important high rate:08_43_49  08_44_10 all:08_21_58 ekf:08_49_36)
ulgFileName = '05_24_43'; % the ulog file name. load 06_28_53 06_45_38  06_54_35 wind 11_39_14  15_57_58 wp:17_09_30, endurance:17_14_37 (18.9) 06_26_32 (17.8) 06_50_20 (19.6) (19.5) 07_10_30 (18.6)
tmp=[ ulgFileName '.mat'];
% exist tmp var05_31_54.ulg
if exist(tmp,"file")
    load(ulgFileName,'log');
else
    if ismac
        % on macOS, run " which ulog2csv " on terminal to get it.
        command = ['!/Users/mch/opt/anaconda3/bin/ulog2csv ' ulgFileName '.ulg']; % /usr/local/bin/ is the path of ulog2csv, 
    else
        % on windows and linux just make sure you have installed pyulog
        command = ['!ulog2csv ' ulgFileName '.ulg']; % have installed ulog2csv,
    end

	eval(command);
    log.data = csv_topics_to_d(ulgFileName);
    log.FileName = ulgFileName;
    log.version = 1.0;
    log.params = '';
    log.messages = '';
    log.info = '';
    %run add_fields_in_preprocessing.m
    save(ulgFileName,'log')
    delete(['*' ulgFileName '*.csv'])
end
%%



if(isfield(log.data, 'allocation_value_0'))
    allocation_value=log.data.allocation_value_0{:,:};

end 
if(isfield(log.data, 'actuator_outputs_value_0'))
    actuator_outputs_value=log.data.actuator_outputs_value_0{:,:};

end 



if(isfield(log.data, 'actuator_outputs_value_0') && isfield(log.data, 'allocation_value_0'))
    % plot((allocation_value(:,1))*1e-6, allocation_value(:,12),'r-','LineWidth',1);hold on;
    plot((allocation_value(:,1))*1e-6, allocation_value(:,16),'k--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    plot((actuator_outputs_value(:,1))*1e-6, actuator_outputs_value(:,2),'b--','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    title('执行器动态');
    xlabel({'时间(s)'});
    ylabel('响应')
    legend('cmd','response','add filter for the same');
end 


