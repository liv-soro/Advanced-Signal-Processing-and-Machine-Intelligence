clear;
%% 1.1 Properties of Power Spectral Density (PSD)
% instantiate a pulse
pulse = zeros(3000, 1);
pulse(1500) = 2;
% ACF of pulse
acf_pulse = xcorr(pulse, 'biased');
% PSD estimates
psd1_pulse = abs(fftshift(fft(acf_pulse)));
psd2_pulse = (abs(fftshift(fft(pulse))).^2)/length(pulse);
f1 = -1000/2:1000/length(acf_pulse):1000/2 - 1000/length(acf_pulse);
f11 = -1000/2:1000/length(psd2_pulse):1000/2 - 1000/length(psd2_pulse);

% instantiate WGN
wgn = randn(3000,1);
% ACF of WGN
acf_wgn = xcorr(wgn, 'biased');
% PSD estimates
psd1_wgn = abs(fftshift(fft(acf_wgn)));
psd2_wgn = (abs(fftshift(fft(wgn))).^2)/length(wgn);
f2 = -1000/2:1000/length(acf_wgn):1000/2 - 1000/length(acf_wgn);
f22 = -1000/2:1000/length(psd2_wgn):1000/2 - 1000/length(psd2_wgn);

% instantiate sinewave
sine = sin(2*pi*0.5*(0:2999)./1000);
% ACF of sinwave
acf_sine = xcorr(sine, 'biased');
% PSD estimates
psd1_sine = abs(fftshift(fft(acf_sine)));
psd2_sine = (abs(fftshift(fft(sine))).^2)/length(sine);
f3 = -1000/2:1000/length(acf_sine):1000/2 - 1000/length(acf_sine);
f33 = -1000/2:1000/length(psd2_sine):1000/2 - 1000/length(psd2_sine);

figure(1)
subplot(3,3,1)
plot(pulse, 'LineWidth',1.3)
title('Pulse', 'FontSize', 14);
xlabel('time (ms)')
ylabel('x(t)')

subplot(3,3,2)
plot(wgn, 'LineWidth',1.5)
title('WGN', 'FontSize', 12);
xlabel('time (ms)')
ylabel('x(t)')

subplot(3,3,3)
plot(sine, 'LineWidth',1.5)
title('Sinewave', 'FontSize', 12);
xlabel('time (ms)')
ylabel('x(t)')

subplot(3,3,4)
plot(f1, acf_pulse, 'LineWidth', 1.5)
title('Pulse ACF', 'FontSize', 12);
xlabel('time lag (k)')
ylabel('ACF')

subplot(3,3,5)
plot(f2, acf_wgn, 'LineWidth', 1.5)
title('WGN ACF', 'FontSize', 12);
xlabel('time lag (k)')
ylabel('ACF')

subplot(3,3,6)
plot(f3, acf_sine, 'LineWidth', 1.5)
title('Sinewave ACF', 'FontSize', 12);
xlabel('time lag (k)')
ylabel('ACF')

subplot(3,3,7)
plot(f1, psd1_pulse);
hold on;
plot(f11, psd2_pulse, '--');
title('Pulse PSD', 'FontSize', 12);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Definition 1','Definition 2');
ytickformat('%.2f')

subplot(3,3,8)
plot(f2, psd1_wgn);
hold on;
plot(f22, psd2_wgn, '--');
title('WGN PSD', 'FontSize', 12);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Definition 1','Definition 2');

subplot(3,3,9)
plot(f3, psd1_sine);
hold on;
plot(f33, psd2_sine, '--');
title('Sinewave PSD', 'FontSize', 12);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Definition 1','Definition 2');
zoom on