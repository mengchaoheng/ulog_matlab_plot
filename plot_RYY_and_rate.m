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

d2r=pi/180;
r2d=180/pi;
%%
ulgFileName = 'log_6_2023-3-6-23-17-18'; % the ulog file name 
tmp=[ulgFileName '.mat'];
% exist tmp var
if exist(tmp,"file")
    load(ulgFileName,'log');
else

    command = ['!/usr/local/bin/ulog2csv ' ulgFileName '.ulg']; % /usr/local/bin/ is the path of ulog2csv, 

    % on macOS, run " which ulog2csv " on terminal to get it.
    % on windows and linux just make sure you have installed pyulog

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
q_0=vehicle_attitude(:,3);
q_1=vehicle_attitude(:,4);
q_2=vehicle_attitude(:,5);
q_3=vehicle_attitude(:,6);
Roll=quat_to_roll(q_0,q_1,q_2,q_3);
Pitch=quat_to_pitch(q_0,q_1,q_2,q_3);
Yaw=quat_to_yaw(q_0,q_1,q_2,q_3);


%% 
figure,
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,3)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Roll rate');
%% 
figure,
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,4)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Pitch rate');
%% 
figure,
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,5)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Yaw rate');

%% and maybe more figure, all in the variable "log.data"
figure,
plot((vehicle_attitude(:,1))*1e-6, Roll*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Roll');%legend('boxoff');

%% and maybe more figure, all in the variable "log.data"
figure,
plot((vehicle_attitude(:,1))*1e-6, Pitch*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Pitch');%legend('boxoff');
%% and maybe more figure, all in the variable "log.data"
figure,
plot((vehicle_attitude(:,1))*1e-6, Yaw*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Yaw');%legend('boxoff');


