% clear all;
close all;
clc;
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
Width  =15;    % (inch)
Height =10;    % 
d2r=pi/180;
r2d=180/pi;
%%
ulgFileName = '02_04_34'; % the ulog file name 
tmp=[ulgFileName '.mat'];
% exist tmp var
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
    vehicle_angular_velocity=log.data.vehicle_angular_velocity_0{:,1:5};
    vehicle_rates_setpoint=log.data.vehicle_rates_setpoint_0{:,:};

    vehicle_attitude=log.data.vehicle_attitude_0{:,:};
    vehicle_attitude_setpoint=log.data.vehicle_attitude_setpoint_0{:,:};
    q_0=vehicle_attitude(:,3);
    q_1=vehicle_attitude(:,4);
    q_2=vehicle_attitude(:,5);
    q_3=vehicle_attitude(:,6);
    Roll=quat_to_roll([q_0 q_1 q_2 q_3]);
    Pitch=quat_to_pitch([q_0 q_1 q_2 q_3]);
    Yaw=quat_to_yaw([q_0 q_1 q_2 q_3]);
    q_d_0=vehicle_attitude_setpoint(:,6);
    q_d_1=vehicle_attitude_setpoint(:,7);
    q_d_2=vehicle_attitude_setpoint(:,8);
    q_d_3=vehicle_attitude_setpoint(:,9);
    Roll_d=quat_to_roll([q_d_0 q_d_1 q_d_2 q_d_3]);
    Pitch_d=quat_to_pitch([q_d_0 q_d_1 q_d_2 q_d_3]);
    Yaw_d=quat_to_yaw([q_d_0 q_d_1 q_d_2 q_d_3]);


    time=vehicle_attitude_setpoint(:,1)*1e-6;
    N= length(time)-1;
    dt=zeros(N,1);
    for i=1:N
        dt(i)=time(i+1)-time(i);
    end


fig1=figure;

plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,3)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,4)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
h=legend('vehicle rates setpoint','vehicle angular velocity');%legend('boxoff');
set(h,'NumColumns',1,'location','northwest');%northwest
set(fig1.CurrentAxes, 'FontSize', 8,'FontName','Times New Roman','LabelFontSizeMultiplier', 1,'TitleFontSizeMultiplier',1,'LineWidth',0.5)
% fig1.CurrentAxes.YAxis.Exponent = -1;
% fig1.CurrentAxes.XTick = [0 1 2 3 4];
% fig1.CurrentAxes.YTick = [-20 0 20 40];
% PlotToFileColorPDF(fig1,'Fig_vehicle_rates.pdf',Width,Height);
%% and maybe more figure, all in the variable "log.data"

fig2=figure;

subplot(3,1,1);
plot((vehicle_attitude(:,1))*1e-6-(vehicle_attitude(1,1))*1e-6, Roll(:)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_attitude_setpoint(:,1))*1e-6-(vehicle_attitude_setpoint(1,1))*1e-6, Roll_d(:)*r2d,'r-','LineWidth',1);hold on;
grid on;
xlabel({'Time(s)'});
ylabel('Roll (deg)')
title('Euler Angle Estimates');
%% and maybe more figure, all in the variable "log.data"
subplot(3,1,2);
plot((vehicle_attitude(:,1))*1e-6-(vehicle_attitude(1,1))*1e-6, Pitch(:)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_attitude_setpoint(:,1))*1e-6-(vehicle_attitude_setpoint(1,1))*1e-6, Pitch_d(:)*r2d,'r-','LineWidth',1);hold on;
grid on;
xlabel({'Time(s)'});
ylabel('Pitch (deg)')
%% and maybe more figure, all in the variable "log.data"
subplot(3,1,3);
plot((vehicle_attitude(:,1))*1e-6-(vehicle_attitude(1,1))*1e-6, Yaw(:)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_attitude_setpoint(:,1))*1e-6-(vehicle_attitude_setpoint(1,1))*1e-6, Yaw_d(:)*r2d,'r-','LineWidth',1);hold on;

