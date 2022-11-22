% 1.3.b)c) Correlation Estimation
clear;

% Generate the PSD estimate (biased) of several realisations of a random
% process

t = [0:0.01:9.99];
reps = 100;
noise = randn(reps, length(t));
x1 = repmat((sin(2*pi*20*t)+sin(2*pi*5*t)), reps, 1) + noise;
L = 2*length(t) - 1;
normf = linspace(-1, 1, L);
P_all = zeros(reps, L);

figure();
for i = 1:reps
    acf = xcorr(x1(i,:), 'biased');
    P = abs(fftshift(fft(acf)));
    
    % Plotting with linear power
    subplot(2, 2, 1); hold on; 
    one = plot(normf, P, 'Color',[0.9290, 0.6940, 0.1250] ); xlim([0, 1]); hold off;
    
    % Plotting with power in dBs
    subplot(2, 2, 3); hold on;
    one_dB = plot(normf, 10*log(P), 'Color', [0.9290, 0.6940, 0.1250]);  xlim([0, 1]); hold off;
    
    % Store PSD estimates
    P_all(i, :) = P;
end
subplot(2, 2, 1); hold on; 
m = plot(normf, mean(P_all, 1), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.2);
set(gca, 'FontSize', 14); xlabel('Normalised Frequency'); ylabel('Power')
title({'PSD Estimates of Sinusoid', 'with Random Noise'}); legend([one, m], {'Realisations', 'Mean'})

subplot(2, 2, 3); hold on; 
m = plot(normf, 10*log10(mean(P_all, 1)), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.2);
set(gca, 'FontSize', 14); xlabel('Normalised Frequency'); ylabel('Power')
title({'PSD Estimates of Sinusoid', 'with Random Noise (dB)'}); legend([one_dB, m], {'Realisations', 'Mean'})

subplot(2, 2, 2); plot(normf, std(P_all, 1), 'LineWidth', 1.2); xlim([0, 1]);
set(gca, 'FontSize', 14); xlabel('Normalised Frequency'); ylabel('Power (dB)')
title('Standard Deviation of PSD Estimate'); 

subplot(2, 2, 4); plot(normf, 10*log10(std(P_all, 1)), 'LineWidth', 1.2); xlim([0, 1]);
set(gca, 'FontSize', 14); xlabel('Normalised Frequency'); ylabel('Power (dB)')
title('Standard Deviation of PSD Estimate (dBs)'); 