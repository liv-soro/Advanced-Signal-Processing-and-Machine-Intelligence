%% 3.2 Adaptive AR Model Based Time-Frequency Estimation
clear; close all; clc;

%% a) Frequency estimation of non-stationary signal with AR(1)

fs = 1000; % sampling frequency;
nSamp = 1500; % length of signal
var_wgn = 0.05; % standard deviation of noise

% Generating the signal
f = [100 * ones(1, 500), ...
     100 + ([501:1000]-500)/2 , ...
     100 + (([1001:1500]-1000)/25).^2]; % frequencies array
phi = cumsum(f); % integral (discrete)
eta = wgn(1, nSamp, pow2db(var_wgn), 'complex'); % complex wgn
y = exp(1i * ((2*pi)/fs * phi)) + eta; % FM signal

% Estimating the order with aryule and the Power Spectrum with freqz
a = aryule(y, 1);
[h, w] = freqz(1, a, nSamp, fs);
P = abs(h).^2;

figure(); hold on; %subplot(2, 2, 1); 
plot(w, 10*log10(P), 'LineWidth', 1.5)
order = 1;

for i = 1:3
    a = aryule(y((i-1)*500 + 1: i*500), order);
    [h, w] = freqz(1, a, nSamp/3, fs);
    P = abs(h).^2;
    %subplot(2, 2, i+1); 
    plot(w, 10*log10(P), 'LineWidth', 1.1)%'Color', [0.8500, 0.3250, 0.0980]
    xlim([0, fs/2])
end
plot(mean(f)*ones(1,1500), linspace(-10,30,1500), '--', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410])

legend('FM Signal', 'Segment 1', 'Segment 2', 'Segment 3', 'Mean frequency')
xlim([0, fs/2])
title('Power Spectrum using AR(1) model', 'FontSize', 13)
xlabel('Frequency (Hz)', 'FontSize', 13); ylabel('Power Spectral Density (dB)', 'FontSize', 13)


%% b) Estimate AR coefficients with CLMS
x = delay(y,1); % input of CLMS = delayed version by of signal (delay = 1)
stepsize = [1, 0.1, 0.01, 0.001];

figure()
for i = 1:length(stepsize)
    [a, ~, ~] = CLMS(x, y, stepsize(i), 1);
    H = zeros(nSamp, nSamp);
    for n = 1:nSamp
        % Run complex-valued LMS algorithm to estimate AR coefficient ahat_1(n)
        [h, w] = freqz(1 , [1; -conj(a(n))], nSamp, fs); % Compute power spectrum
        H(:, n) = abs(h).^2; % Store it in a matrix 
    end
    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    % Plot time-frequency diagram
    subplot(4, 1, i); 
    mesh(1:nSamp, w, H)
    view(2)
    xlabel('Time Step', 'Fontsize', 13); ylabel('Frequency (Hz)', 'Fontsize', 13);
    cb = colorbar;
    ylabel(cb, 'PSD (dB)', 'Fontsize', 13);
    title("Spectral Estimate with Adaptive AR Coefficient (CLMS with \mu =  " + stepsize(i) + ")", 'Fontsize', 13)
end



