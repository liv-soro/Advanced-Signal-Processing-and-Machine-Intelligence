% 1.5 Real World Signals: Respiratory Sinus Arrhythmia from RR-intervals

% Load Data
load('RAW-ECG.mat')
load('RRI-DATA.mat')

% Find sampling frequency of fRRIs
 Trial_1 = data(770:2.5e5);
 [xRRI1,fsRRI1]=ECG_to_RRI(Trial_1, fs);

%% 1.5.a)
% Preprocessing
xRRI1_mean = mean(xRRI1);
xRRI1_detrend = detrend(xRRI1-xRRI1_mean);
xRRI2_mean = mean(xRRI2);
xRRI2_detrend = detrend(xRRI2-xRRI2_mean);
xRRI3_mean = mean(xRRI3);
xRRI3_detrend = detrend(xRRI3-xRRI3_mean);
xRRI = {xRRI1_detrend, xRRI2_detrend, xRRI3_detrend};
%fsRRI = fsRRI1;
fsRRI = 4;

% Parameters
noverlap = 0;
w_std = cell(1, length(xRRI));
w_avg = cell(3, length(xRRI));

% PSD Estimates
PSD_std = cell(1, length(xRRI));
PSD_avg = cell(3, length(xRRI));

figure();
% Averaged Periodogram (Bartlett's method)
for i_RRI = 1:length(xRRI)
    n = length(xRRI{i_RRI});
    wind_n = [n, 50*fsRRI, 150*fsRRI];
    for i_wind = 1:length(wind_n)
        [PSD, w] = pwelch(xRRI{1, i_RRI}, wind_n(1, i_wind), noverlap, length(xRRI{1, i_RRI}), fsRRI);
        subplot(3, 1, i_RRI)
        plot(w, 10*log10(PSD), 'LineWidth', 1.3); hold on
    end
    legend('Standard', 'Window = 50s', 'Window = 150s')
    title(sprintf('PSD estimates of RRI %d', i_RRI))
    ylabel('Power (dBs)'); xlabel('Frequency (Hz)')
    set(gca, 'Fontsize', 13);
end

%% Confirming frequency peaks with MUSIC

%% 1.5.c)

% Determining model order with pyulear
order = 1:10;
figure();
for i_RRI = 1:3
    subplot(3,1,i_RRI)
    for o = 1: length(order)
        [PSD, w] = pyulear(xRRI{i_RRI}, order(o), 2048, fsRRI);
        plot(w, 10*log10(PSD)); hold on
    end
    legend('1', '2', '3', '4','5', '6', '7','8', '9', '10')
end

%% Determining model order with aryule
order = 1:10;
figure();
for i_RRI = 1:3
    subplot(3,1,i_RRI)
    n = length(xRRI{i_RRI});
    for o = 1: length(order)
        [ar_coeffs, var] = aryule(xRRI{i_RRI}, order(o));
        [h, w] = freqz(sqrt(var), ar_coeffs, n, fsRRI);
        plot(w, abs(h).^2); hold on
    end
    legend('1', '2', '3', '4','5', '6', '7','8', '9', '10')
end
%% Plotting With Orders found
figure();
order = [2, 9, 3]; %[3, 9, 4] % for trial 1 say that you can't find model order 
%(unconstrained breathing)
for i_RRI = 1:3
    [PSD, w] = pyulear(xRRI{i_RRI}, order(i_RRI), 2048, fsRRI);
    plot(w, 10*log10(PSD), 'LineWidth', 1.4); hold on
end
legend('Trial 1, Order 2', 'Trial 2, Order 9', 'Trial 3, Order 3')
set(gca, 'FontSize', 13)
xlabel('Frequency (Hz)'); ylabel('Power (dBs)')
title('Autoregressive PSD estimate of RRI data')



