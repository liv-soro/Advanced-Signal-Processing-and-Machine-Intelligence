% 1.2.a) Periodogram-based methods applied to real-world data
load sunspot.dat

% Preprocessing: remove mean and detrend
sunspot_mean = mean(sunspot(:,2));
sunspot_detrend = detrend(sunspot(:,2)-sunspot_mean);
figure(); 
plot(sunspot(:,1), sunspot(:,2), 'LineWidth', 1.4); hold on
plot(sunspot(:,1), sunspot_detrend,'LineWidth', 1.4)

% Preprocessing: log
sunspot_log = log(sunspot(:,2) + eps); % add eps to prevent taking the log of 0
sunspot_log_mean = sunspot_log - mean(sunspot_log);

plot(sunspot(:,1), sunspot_log_mean, 'LineWidth', 1.4);
xlim([1700, 1987])
ylim([-60,200])
legend('Original','Mean-Detrend','Log-Mean') %change legend
title('Preprocessing', 'FontSize', 14) %change title (maybe not needed)
xlabel('Years', 'FontSize', 14)
ylabel('Signal Values', 'FontSize', 14)


% Peridodgram-based spectral estimation
[psd_raw, w] = periodogram(sunspot(:,2), hamming(length(sunspot(:,2))));
figure();
plot(w/pi, log10(psd_raw).*10, 'LineWidth', 1.3); hold on
psd_detrend = periodogram(sunspot_detrend, hamming(length(sunspot(:,2))));
plot(w/pi, log10(psd_detrend).*10, '--', 'LineWidth', 1.3)
psd_log = periodogram(sunspot_log_mean, hamming(length(sunspot(:,2))));
plot(w/pi, log10(psd_log).*10, 'LineWidth', 1.3)
xlim([0, 1])
ylim([-30, 60])
legend('Original','Mean & Detrend','Log-Mean')
title('PSD Estimates for Sunspot Time Series', 'FontSize', 14) 
xlabel('Normalised Frequency (\pi rad/sample)', 'FontSize', 14)
ylabel('Power (dBs)', 'FontSize', 14)
% log is data compression so the peaks seen in the log are actually the
% significant peaks

