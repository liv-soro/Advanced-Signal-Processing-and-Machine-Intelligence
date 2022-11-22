% Properties of the Power Spectral Desnity (PSD)

Fs = 1000; % sampling frequency
T = 1/Fs; % sampling period
L = 2000; % length of signal
t = (0:L-1)./Fs; % time vector

% Sinusoid
x1 = sin(2*pi*0.5*t); % sine wave
rxx1 = xcorr(x1, 'biased'); % autocorrelation
P_rxx1 = abs(fftshift(fft(rxx1))); % def 1 %why abs?
f1 = -Fs/2 : Fs/length(rxx1) : Fs/2 - Fs/length(rxx1); % frequency axis

P_est = (abs(fftshift(fft(x1))).^2)/L; % def 2
f2 = -Fs/2 : Fs/length(P_est) : Fs/2 - Fs/length(P_est); % frequency axis

% Pulse
x2 = zeros(L, 1);
x2(L/2) = 10; % pulse
rxx2 = xcorr(x2, 'biased'); % autocorrelation
P_rxx2 = abs(fftshift(fft(rxx2))); % def 1
f3 = -Fs/2 : Fs/length(rxx2) : Fs/2 - Fs/length(rxx2); % frequency axis

P_est2 = (abs(fftshift(fft(x2))).^2)/L; % def 2
f4 = -Fs/2 : Fs/length(P_est2) : Fs/2 - Fs/length(P_est2); % frequency axis



% Plots
figure();
subplot(2, 3, 1); plot(t, x2, 'LineWidth', 1.3);
set(gca, 'FontSize', 14); xlabel('t(s)'); ylabel('Amplitude'); title('Pulse');

subplot(2, 3, 2); plot([-L+1 : L-1], rxx2, 'LineWidth', 1.3)
set(gca, 'FontSize', 14); xlabel('Time Lag (k)'); ylabel('ACF'); title('Fast Decaying ACF');

subplot(2, 3, 3); plot(f3, P_rxx2, 'LineWidth', 1.3); hold on; plot(f4, P_est2, '--', 'LineWidth', 1.3)
set(gca, 'FontSize', 14); legend('Def 1', 'Def 2'); %xlim([-10,10]); ylim([0,0.1])
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('PSD Estimates');

subplot(2, 3, 4); plot(t, x1, 'LineWidth', 1.3);
set(gca, 'FontSize', 14); xlabel('Time (s)'); ylabel('Amplitude'); title('Sinusoidal Signal');

subplot(2, 3, 5); plot([-L+1 : L-1], rxx1, 'LineWidth', 1.3)
set(gca, 'FontSize', 14); xlabel('Time Lag (k)'); ylabel('ACF'); title('Slow Decaying ACF');

subplot(2, 3, 6); plot(f1, P_rxx1, 'LineWidth', 1.3); hold on; plot(f2, P_est, '--', 'LineWidth', 1.3)
set(gca, 'FontSize', 14); legend('Def 1', 'Def 2');% xlim([-10,10])
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('PSD Estimates');



%%
figure();
subplot(2, 3, 1); plot(pulse, 'LineWidth', 1.3);
set(gca, 'FontSize', 14); xlabel('t(s)'); ylabel('Amplitude'); title('Pulse');

subplot(2, 3, 2); plot(f1, acf_pulse, 'LineWidth', 1.3)
set(gca, 'FontSize', 14); xlabel('Time Lag (k)'); ylabel('ACF'); title('Fast Decaying ACF');

subplot(2, 3, 3); plot(f1, psd1_pulse, 'LineWidth', 1.3); hold on; plot(f11, psd2_pulse, '--', 'LineWidth', 1.3)
set(gca, 'FontSize', 14); legend('Def 1', 'Def 2'); 
%xlim([-10,10]); 
ylim([0,0.01])
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('PSD Estimates');

subplot(2, 3, 4); plot(sine, 'LineWidth', 1.3);
set(gca, 'FontSize', 14); xlabel('Time (s)'); ylabel('Amplitude'); title('Sinusoidal Signal');

subplot(2, 3, 5); plot(f3, acf_sine, 'LineWidth', 1.3)
set(gca, 'FontSize', 14); xlabel('Time Lag (k)'); ylabel('ACF'); title('Slow Decaying ACF');

subplot(2, 3, 6); plot(f3, psd1_sine, 'LineWidth', 1.3); hold on; plot(f33, psd2_sine, '--', 'LineWidth', 1.3)
set(gca, 'FontSize', 14); legend('Def 1', 'Def 2');% xlim([-10,10])
xlabel('Frequency (Hz)'); ylabel('Magnitude'); title('PSD Estimates');


%---------------------------------------%
% % WGN
% x2 = wgn(L, 1, 1); % white gaussian noise
% rxx2 = xcorr(x2, 'biased'); % autocorrelation
% P_rxx2 = abs(fftshift(fft(rxx2))); % def 1
% f3 = -Fs/2 : Fs/length(rxx2) : Fs/2 - Fs/length(rxx2); % frequency axis
% 
% P_est2 = (abs(fftshift(fft(x2))).^2)/L; % def 2
% f4 = -Fs/2 : Fs/length(P_est2) : Fs/2 - Fs/length(P_est2); % frequency axis

