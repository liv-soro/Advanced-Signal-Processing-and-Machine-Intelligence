% Adaptive Noise Cancellation
clc; clear; close all

%% Generate noise corrupted signal
t = 1:1000;
x = sin(0.01*pi*t); x = x'; % clean signal

eta = filter([1, 0, 0.5], [1], randn(1500, 1)); % coloured noise
eta = eta(501:end); % remove transient

s = x + eta; % input signal (noise corrupted signal)
% figure(); plot(t, s)

L = length(s);

%% a) ALE - Justifying choice for minimum delta
% Confirm theoretical result of minimum value for the delay
s_R = xcorr(s, 'biased');
figure(); 
subplot(1, 2, 1); hold on
plot([-(nSamp-1):(nSamp-1)], s_R, 'LineWidth', 1.3) % check that after lag = 2, the ACF is dominated by the sinusoid
plot(3*ones(1,1500), linspace(-0.5,2,1500), 'LineWidth', 1.5, 'Color', [0.4660, 0.6740, 0.1880])
xlim([-100,100])
xlabel('Time Lag (k)'); ylabel('ACF')
title('ACF of s(n)'); legend('ACF', '\Delta_{min}')
set(gca, 'FontSize', 13)


% MSPE vs Filter Length & Delay
M = 5:5:20; % filter lenghts
delta = 1:25;
stepsize = 0.01;
MSPE = zeros(length(M), length(delta));
leg = cell(1,length(M));

%figure()
subplot(1, 2, 2); 
hold on
for m = 1:length(M)
    for d = 1:length(delta)
        [~, xhat, ~] = LMS_ALE(s, stepsize, delta(d), M(m));
        MSPE(m, d) = (1/L) * dot((x-xhat), (x-xhat));
    end
    plot(delta, 10*log10(MSPE(m, :)), 'LineWidth', 1.4);
    leg{m} = ("M = " + num2str(M(m)));
end
plot(3*ones(1,1500), linspace(-6,-2,1500),  'LineWidth', 1.5)
leg{m+1} = '\Delta_{min}';
legend(leg, 'FontSize', 13)
xlabel('Delay (\Delta)', 'FontSize', 13); ylabel('MSPE (dB)', 'FontSize', 13);
title('MSPE vs Filter Length & Delay', 'FontSize', 13);


%% Comparing the signals
figure()
[~, xhat3, ~] =  LMS_ALE(s, stepsize, 3, 5);
subplot(1, 2, 2); hold on;
plot(1:length(s), s', 'LineWidth', 1)
plot(1:length(xhat), xhat3', 'LineWidth', 1)
plot(1:length(x), x', 'LineWidth', 1.3)
xlabel('Time Step (n)', 'FontSize', 13); ylabel('Amplitude (AU)', 'FontSize', 13);
legend('Noise-corrupted Signal', 'Denoised Signal (ALE)', 'Clean Signal')
title('ALE Signal Denoising', 'FontSize', 13);
%set(gca, 'FontSize', 13);
ylim([-5, 5])

%% b) ALE - Effect of filter order M and delay Delta on ALE

% MSPE vs Delay
M = 5; % filter lenghts
delta = 3:25;
stepsize = 0.01;

MSPE= zeros(length(delta), 1);

for d = 1:length(delta)
    [~, xhat, ~] = LMS_ALE(s, stepsize, delta(d), M);
    MSPE(d, 1) = (1/L) * dot(x-xhat, x-xhat);
end

figure(); subplot(1, 2, 1); hold on
plot(delta, 10*log10(MSPE), 'LineWidth', 2);
xlabel('Delay (\Delta)'); ylabel('MSPE (dB)');
title('MSPE vs Delay (M = 5)');
set(gca, 'FontSize', 13);

% MSPE vs Filter Length
M = 5:5:20; % filter lenghts
delta = 3;
MSPE= zeros(length(M), 1);

for m = 1:length(M)
    [~, xhat, ~] = LMS_ALE(s, stepsize, delta, M(m));
    MSPE(m, 1) = (1/L) * dot(x-xhat, x-xhat);
end

subplot(1, 2, 2); hold on
plot(M, 10*log10(MSPE), 'LineWidth', 1.4, 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Filter Length (M)'); ylabel('MSPE (dB)');
title('MSPE vs Filter Length (\Delta = 3)');
set(gca, 'FontSize', 13);

%% c) Comparing ANC and ALE
clear; close all; clc;
% run code to generate the signals
%%
stepsize = 0.01;
sec_noise = eta + 0.5 * delay(eta, 1)+ 0.2; % Secondary noise is correlated in some unknown way to the primary noise

% MSPE vs Filter Length (ANC)
M = 1:20; % filter lenghts
MSPE = zeros(length(M), 1);

for m = 1:length(M)
    [~, xhat, ~] = LMS_ANC(s, sec_noise, stepsize, M(m));
    MSPE(m, 1) = (1/L) * dot(x-xhat, x-xhat);
end

figure();
plot(M, 10*log10(MSPE), 'LineWidth', 1.4, 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Filter Length (M)'); ylabel('MSPE (dB)');
title('MSPE vs Filter Length ANC');
set(gca, 'FontSize', 13);

[~, xhat_ALE, ~] = LMS_ALE(s, stepsize, 3, 5);
[~, xhat_ANC, ~] = LMS_ANC(s, sec_noise, stepsize, 4);

% Comparing the signals
figure(); subplot(1, 2, 1); hold on
plot(1:length(s), s')
plot(1:length(xhat), xhat_ALE')
plot(1:length(x), x', 'LineWidth', 1.5)
xlabel('Time Step (n)', 'FontSize', 13); ylabel('Amplitude (AU)', 'FontSize', 13);
legend('Noise-corrupted Signal', 'Denoised Signal (ALE)', 'Clean Signal')
title('ALE Signal Denoising', 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(s), s')
pause(1)
plot(1:length(xhat), xhat_ANC')
pause(1)
plot(1:length(x), x', 'LineWidth', 1.5)

xlabel('Time Step (n)', 'FontSize', 13); ylabel('Amplitude (AU)', 'FontSize', 13);
legend('Noise-corrupted Signal', 'Denoised Signal (ALE)', 'Clean Signal')
title('ANC Signal Denoising', 'FontSize', 13);


MSPE_ALE = 10*log10(movmean(1/L * (x-xhat_ALE).^2, 20));
MSPE_ANC = 10*log10(movmean(1/L * (x-xhat_ANC).^2, 20));

figure(); hold on;
plot(1:length(MSPE_ALE), MSPE_ALE, 'LineWidth', 1.3);
plot(1:length(MSPE_ANC), MSPE_ANC, 'LineWidth', 1.3);
xlabel('Time Step (n)', 'FontSize', 13); ylabel('MSPE (dB)', 'FontSize', 13); 
title('MSPE: ALE vs ANC', 'FontSize', 13);
legend('ALE', 'ANC');


%% d) ANC for 50Hz mains noise
clear, close all; clc
%%
%load('EEG_Data_Assignment2.mat')
POz_mean = mean(POz);
data = detrend(POz-POz_mean);

t = [1:length(POz)]';
sec_noise = sin(2*pi*(50/fs)*t) + randn(length(data),1); % reference input

stepsize = [0.1, 0.01, 0.001];
M = 7:11;
nfft = 2^14;

figure();
count = 0;
for i = 1:length(stepsize)
    for j = 1:length(M)
        count = count + 1;
        [~, xhat, ~] = LMS_ANC(data, sec_noise, stepsize(i), M(j));
        subplot(length(stepsize), length(M), count);
        spectrogram(xhat, rectwin(3540), round(0.333*(3540)), nfft, fs, 'yaxis')
        ylim([0, 55])
    end
end

%%
figure();
subplot(1, 2, 1);
spectrogram(data, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis');
ylim([0, 55])
title('Spetrogram of POz Data')
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.01, 40); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title('Spetrogram of Cleaned POz Data (M = 40, \mu = 0.01)')
set(gca, 'FontSize', 12);

%%
figure();
subplot(1, 3, 1);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.01, 20); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title({'Spectrogram of Cleaned POz Data', '(M = 20, \mu = 0.01)'})
set(gca, 'FontSize', 13);

subplot(1, 3, 2);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.01, 40); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title({'Spectrogram of Cleaned POz Data', '(M = 40, \mu = 0.01)'})
set(gca, 'FontSize', 13);

subplot(1, 3, 3);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.01, 60); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title({'Spectrogram of Cleaned POz Data', '(M = 60, \mu = 0.01)'})
set(gca, 'FontSize', 13);

%%
figure();
subplot(1, 3, 1);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.001, 40); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title({'Spectrogram of Cleaned POz Data', '(M = 40, \mu = 0.001)'})
set(gca, 'FontSize', 13);

subplot(1, 3, 2);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.01, 40); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title({'Spectrogram of Cleaned POz Data', '(M = 40, \mu = 0.01)'})
set(gca, 'FontSize', 13);

subplot(1, 3, 3);
[~, xhat, ~] = LMS_ANC(POz, sec_noise, 0.03, 40); % 11; 35; 40
spectrogram(xhat, hamming(3540), round(0.333*(3540)), nfft, fs, 'yaxis') %nfft
ylim([0, 55])
title({'Spectrogram of Cleaned POz Data', '(M = 40, \mu = 0.03)'})
set(gca, 'FontSize', 13);