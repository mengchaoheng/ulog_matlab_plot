clear all;
close all;
clc;
addpath(genpath(pwd));

X =load('data1.txt');

L = size(X,1);   % Length of signal
L = size(X,1);   % Length of signal
T=0.0013; % Sampling period  (second)
t = (0:L-1)*T;        % Time vector

M_PI_F=3.14159265;
sample_freq=1000;
notch_freq=100;
bandwidth=100;

alpha = tan(M_PI_F * bandwidth / sample_freq);
beta = -cos(2 * M_PI_F * notch_freq / sample_freq);
a0_inv = 1 / (alpha + 1);

b0 = a0_inv;
b1 = 2 * beta * a0_inv;
b2 = a0_inv;
a1 = b1;
a2 = (1 - alpha) * a0_inv;

delay_element_1 = X(1,1);
delay_element_2 = delay_element_1;
delay_element_output_2 = -0.035260253 * (b0 + b1 + b2) / (1 + a1 + a2);
delay_element_output_1 = delay_element_output_2;

X_Num=size(X(:,1));
output = X;
for i=1:X_Num
    sample=X(i,1);
    output(i,1) = b0 * sample + b1 * delay_element_1 + b2 * delay_element_2 - a1 * delay_element_output_1 - a2 * delay_element_output_2;
    delay_element_2 = delay_element_1;
    delay_element_1 = sample;
    delay_element_output_2 = delay_element_output_1;
    delay_element_output_1 = output(i,1);
end
figure,
    plot(t,output(:,1));
    hold on;
    plot(t,X(:,1));

title('Signal')
xlabel('t (seconds)')
ylabel('X(t)')
legend('X fil','X','Location','northeastoutside','box','off');

Y = fft(X);
Y_fil = fft(output);
P2 = abs(Y/L);
P1 = P2(1:L/2+1,1);
P1(2:end-1,1) = 2*P1(2:end-1,1);

P2_fil = abs(Y_fil/L);
P1_fil = P2_fil(1:L/2+1,1);
P1_fil(2:end-1,1) = 2*P1_fil(2:end-1,1);


figure,
f = sample_freq*(0:(L/2))/L;
plot(f,P1(:,1));
hold on;
plot(f,P1_fil(:,1))
title('FFT of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('Y fil','Y','Location','northeastoutside','box','off');
