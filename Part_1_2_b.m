% 1.2.b) Periodogram-based methods applied to real-world data

% Load data
EEG_Data = load('EEG_Data_Assignment1');
fs = EEG_Data.fs;
EEG = EEG_Data.POz;

% Preprocessing
EEG_mean = mean(EEG);
EEG_detrend = detrend(EEG-EEG_mean);

% Periodogram
ndft = 5*fs; % 5 DFT samples per Hz
[psd_EEG, w] = periodogram(EEG_detrend, hamming(length(EEG)), ndft, fs);
figure(); subplot(1, 2, 1)
plot(w, log10(psd_EEG).*10, 'LineWidth', 1.3)
xlim([0,60]) % frequency range we are interested in
ylim([-160, -90])
title('Periodogram for EEG data');  set(gca, 'FontSize', 14);
xlabel('Frequency (Hz)')
ylabel('Power (dBs)')

wind = [10, 5, 1];
wind_n = wind.*fs;
noverlap = 0;
subplot(1, 2, 2)
for i = 1:3
    [psd, w] = pwelch(EEG_detrend, wind_n(i), noverlap, length(EEG), fs);
    plot(w, log10(psd).*10, 'LineWidth', 1.3); hold on
end
xlim([0,60]); ylim([-160, -90]);
set(gca, 'FontSize', 14);
leg = legend('10 s','5 s','1 s'); title(leg, 'Window Length');
title('PSD estimates for Windowed EEG data')
xlabel('Frequency (Hz)')
ylabel('Power (dBs)')
%saveas(gcf, '1_2_b(2)', 'epsc')

% Peaks corresponding to SSVEP: peak at 13 and harmonic (26) X = 13
% Compare std approach with 10s windonw: in 10s apporach it is easier to
% distinguish estimated ssvep - why
% Effect of small window size: smaller variance, higher bias

