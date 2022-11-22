%% 3.3 A real time Spectrum Analyser using LMS
clear; close all; clc;

%% c) DFT-CLMS
fs = 1000; % sampling frequency;
nSamp = 1500; % length of signal
var_wgn = 0.05; % standard deviation of noise

% Generating the (time-domain) FM signal
f = [100 * ones(1, 500), ...
     100 + ([501:1000]-500)/2 , ...
     100 + (([1001:1500]-1000)/25).^2]; % frequencies array
phi = cumsum(f); % integral (discrete)
eta = wgn(1, nSamp, pow2db(var_wgn), 'complex'); % complex wgn
y = exp(1i * ((2*pi)/fs * phi)) + eta; % FM signal

%% Generating the Fourier signal
nSamp = 1500;
n = 1:nSamp;
expo = @(n) (exp((1i * 2 * n * (0:nSamp - 1) * pi)/nSamp));
x = expo((1:nSamp)'); x = x/nSamp;

%% DFT-CLMS on FM signal
x_in = delay(y, 1);

figure();
leaks = [0, 0.01, 0.1, 0.5];
for j = 1:length(leaks)
    [h_clms, ~, ~] = CLMS_dft(x, y, 1, leaks(j));
    subplot(length(leaks), 1, j)
    mesh(abs(h_clms).^2)
    view(2)
    xlabel('Time Step', 'Fontsize', 13); ylabel('Frequency (Hz)', 'Fontsize', 13);
    cb = colorbar;
    ylabel(cb, 'PSD (dB)', 'Fontsize', 13);
    xlim([0, 1500]); ylim([0, 1000])
    title("Spectral Estimate with DFT-CLMS (\gamma = " + leaks(j) + ")", 'Fontsize', 13)
end


%% DFT-CLMS on EEG signal
% Load Data
addpath('Data');
load('EEG_Data_Assignment2.mat');

POz_mean = mean(POz);
POz = detrend(POz - POz_mean); % remove mean and detrend

nSamp = 1200;
f = 0:fs-1;
a = 4000; % choose a starting point for the sebment
t = a:a+nSamp-1;
POz = POz(t);

expo = @(n) (exp((1i * 2 * n * (0:nSamp - 1) * pi)/nSamp));
x = expo((1:nSamp)'); x = x/nSamp;

[h_clms, ~, ~] = CLMS_dft(x, POz, 1, 0);

figure(); 
mesh(t, f, abs(h_clms(:, 1:end-1)).^2) % careful on how the time axis appears
view(2)
ylim([0,60])
cb = colorbar;
ylabel(cb, 'PSD (dB)', 'Fontsize', 13);
xlabel('Time step'); ylabel('Frequency (Hz)');
title('Spectral Estimate of EEG POz signal with DFT-CLMS')
set(gca, 'FontSize', 13)


