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
%%
ulgFileName = '09_07_56'; % the ulog file name 
tmp=[ ulgFileName '.mat'];
% exist tmp var
if exist(tmp,"file")
    load(ulgFileName,'log');
else
    if ismac
        % on macOS, run " which ulog2csv " on terminal to get it.
        command = ['!/usr/local/bin/ulog2csv ' ulgFileName '.ulg']; % /usr/local/bin/ is the path of ulog2csv, 
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
%     delete(['*' ulgFileName '*.csv'])
end


vehicle_angular_velocity=log.data.vehicle_angular_velocity_0{:,1:5};
vehicle_rates_setpoint=log.data.vehicle_rates_setpoint_0{:,:};
vehicle_attitude=log.data.vehicle_attitude_0{:,:};
vehicle_attitude_setpoint=log.data.vehicle_attitude_setpoint_0{:,:};
q_0=vehicle_attitude(:,3);
q_1=vehicle_attitude(:,4);
q_2=vehicle_attitude(:,5);
q_3=vehicle_attitude(:,6);
Roll=quat_to_roll(q_0,q_1,q_2,q_3);
Pitch=quat_to_pitch(q_0,q_1,q_2,q_3);
Yaw=quat_to_yaw(q_0,q_1,q_2,q_3);
q_0_setpoint=vehicle_attitude_setpoint(:,6);
q_1_setpoint=vehicle_attitude_setpoint(:,7);
q_2_setpoint=vehicle_attitude_setpoint(:,8);
q_3_setpoint=vehicle_attitude_setpoint(:,9);
Roll_setpoint=quat_to_roll(q_0_setpoint,q_1_setpoint,q_2_setpoint,q_3_setpoint);
Pitch_setpoint=quat_to_pitch(q_0_setpoint,q_1_setpoint,q_2_setpoint,q_3_setpoint);
Yaw_setpoint=quat_to_yaw(q_0_setpoint,q_1_setpoint,q_2_setpoint,q_3_setpoint);

rate_N=size(vehicle_rates_setpoint(:,1));
rate_delta_t=zeros(rate_N-1);
for i=1:rate_N-1
rate_delta_t(i)=(vehicle_rates_setpoint(i+1,1))*1e-6-(vehicle_rates_setpoint(i,1))*1e-6;
end

attitude_N=size(vehicle_attitude_setpoint(:,1));
attitude_delta_t=zeros(attitude_N-1);
for i=1:attitude_N-1
attitude_delta_t(i)=(vehicle_attitude_setpoint(i+1,1))*1e-6-(vehicle_attitude_setpoint(i,1))*1e-6;
end

data_LM=load('all_10.00.50(1).csv')/100;

len =size(data_LM);
Roll_LM=data_LM(:,1);
Pitch_LM=data_LM(:,2);
Yaw_LM=data_LM(:,3);
dt=0.023;
t=0:dt:(len-1)*dt;

start=450;
start2=1;
%% and maybe more figure, all in the variable "log.data"
figure,
plot(t(start:end)-t(start), Roll_LM(start:end)*1+2,'k-','LineWidth',1);hold on;

plot((vehicle_attitude(start2:end,1))*1e-6-(vehicle_attitude(start2,1))*1e-6, Roll(start2:end)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Roll');%legend('boxoff');

%% and maybe more figure, all in the variable "log.data"
figure,
plot(t(start:end)-t(start), Pitch_LM(start:end)*1,'k-','LineWidth',1);hold on;

plot((vehicle_attitude(start2:end,1))*1e-6-(vehicle_attitude(start2,1))*1e-6, Pitch(start2:end)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Pitch');%legend('boxoff');
%% and maybe more figure, all in the variable "log.data"
figure,
plot(t(start:end)-t(start), Yaw_LM(start:end)*1,'k-','LineWidth',1);hold on;

plot((vehicle_attitude(start2:end,1))*1e-6-(vehicle_attitude(start2,1))*1e-6, Yaw(start2:end)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Angular velocity(deg/s)')
% title('Pitch angular rate');
legend('Yaw');%legend('boxoff');