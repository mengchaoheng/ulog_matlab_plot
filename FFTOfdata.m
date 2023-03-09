%% Noisy Signal
% Use Fourier transforms to find the frequency components of a signal buried 
% in noise.
% 假设数据按列存储，每一列是一种数据，L为数据长度，num为存储数据的数量
X=load('all_17.01.22.csv');
%% 
% calculate the sampling time and frequency.

T=0.018; % Sampling period  (second)
Fs=1/T; % Sampling frequency
L = size(X,1);   % Length of signal
t = (0:L-1)*T;        % Time vector
num = size(X,2);   % num of data type, like acc_x,acc_y,acc_z, then the num = 3
le_str = cell(num,1);  % for legen

%% 
% 
%% 
% Plot the noisy signal in the time domain. It is difficult to identify the 
% frequency components by looking at the signal |X(t)|. 

figure,
for i = 1:num
    plot(t,X(:,i));
    le_str{i} = ['X_',num2str(i)];
    hold on;
end

title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (seconds)')
ylabel('X(t)')
legend(le_str,'Location','northeastoutside','box','off');





%% 
% Compute the Fourier transform of the signal. 

Y = fft(X);
%% 
% Compute the two-sided spectrum |P2|. Then compute the single-sided spectrum 
% |P1| based on |P2| and the even-valued signal length |L|.

P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
%% 
% Define the frequency domain |f| and plot the single-sided amplitude spectrum 
% |P1|. The amplitudes are not exactly at 0.7 and 1, as expected, because of the 
% added noise. On average, longer signals produce better frequency approximations.
figure,
f = Fs*(0:(L/2))/L;
for i = 1:num
    plot(f,P1(:,i));
    le_str{i} = ['P1_',num2str(i)];
    hold on;
end
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend(le_str,'Location','northeastoutside','box','off');
