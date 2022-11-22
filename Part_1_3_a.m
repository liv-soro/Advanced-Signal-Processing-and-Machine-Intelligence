% 1.3.a) Correlation Estimation

t = [0:0.01:9.99];
fs = 100;
%% WGN - ACF & Correlograms
x1 = randn(size(t));
%x1 = wgn(length(t), 1, 1);

acf1_unb = xcorr(x1, 'unbiased'); acf1_b = xcorr(x1, 'biased');

P1_unb = abs(fftshift(fft(acf1_unb))); P1_b = abs(fftshift(fft(acf1_b)));

normf = linspace(-1,1, length(acf1_unb)); %dunno about this...
% is it normal that I don't use any fs

%% Noisy sinusoid - ACF and Correlograms
x2 = 5*sin(2*pi*20*t) + x1;
acf2_unb = xcorr(x2, 'unbiased');
acf2_b = xcorr(x2, 'biased');

P2_unb = abs(fftshift(fft(acf2_unb)));
P2_b = abs(fftshift(fft(acf2_b)));

%% AR(4) model - ACF and Correlograms
mdl = arima('Constant', 0, 'AR', {2.21, -2.94, 2.17, -0.96}, 'Variance', 1); % coefficients from lecture notes
x3 = simulate(mdl, length(t));

acf3_unb = xcorr(x3, 'unbiased');
acf3_b = xcorr(x3, 'biased');

P3_unb = abs(fftshift(fft(acf3_unb)));
P3_b = abs(fftshift(fft(acf3_b)));

%% Plots of the ACF
% ACF of WGN
figure();
subplot(1, 3, 1);
plot(1:length(acf1_unb), acf1_unb); hold on
plot(1:length(acf1_b), acf1_b)
title('ACF estimates of wgn', 'FontSize', 13)
xlabel('Correlation Lag', 'FontSize', 13)
ylabel('ACF', 'FontSize', 13)
legend('Unbiased', 'Biased')
xlim([0, length(acf1_b)])

% ACF of Noisy Sinusoid
subplot(1, 3, 2)
plot(1:length(acf2_unb), acf2_unb); hold on
plot(1:length(acf2_b), acf2_b)
title('ACF estimates of a noisy sinusoid', 'FontSize', 13)
xlabel('Correlation Lag', 'FontSize', 13)
ylabel('ACF', 'FontSize', 13)
legend('Unbiased', 'Biased')
xlim([0, length(acf2_b)])

% Filtered WGN
a = [0.5, -0.5];
a = [1, -a];
x4 = filter(1, a, x1);
acf4_unb = xcorr(x4, 'unbiased');
acf4_b = xcorr(x4, 'biased');

subplot(1, 3, 3)
plot(1:length(acf4_unb), acf4_unb); hold on
plot(1:length(acf4_b), acf4_b)
title('ACF estimates of filtered wgn', 'FontSize', 13)
xlabel('Correlation Lag', 'FontSize', 13)
ylabel('ACF', 'FontSize', 13)
legend('Unbiased', 'Biased')
xlim([0, length(acf4_b)])

% the estimates differ more and more with increasing k (lag)
% as the unbiased estimator start having increasing values: 
% -> data available more outdated
% in the biased estimator: you attenuate the bad/old estimates

%see lecture 3, slide 19 (positive semidefinite)
%% Plots of the Correlograms:
figure()
subplot(3, 2, 1); plot(normf, 10*log10(P1_unb)); 
xlim([0, 1]); % plotting only for positive frequencies
xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14, 'ylim', [-40,20]); title('PSD of WGN - Unbiased');

subplot(3, 2, 2); plot(normf, 10*log10(P1_b)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14, 'ylim', [-40,20]); title('PSD of WGN - Biased');

subplot(3, 2, 3); plot(normf, 10*log10(P2_unb)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14, 'ylim', [-20,40]); title('PSD of Noisy Sinusoid - Unbiased');

subplot(3, 2, 4); plot(normf, 10*log10(P2_b)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14, 'ylim', [-20,40]); title('PSD of Noisy Sinusoid - Biased');

subplot(3, 2, 5); plot(normf, 10*log10(P3_unb)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14, 'ylim', [-50,50]); title('PSD of AR(4) Process - Unbiased');

subplot(3, 2, 6); plot(normf, 10*log10(P3_b)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14, 'ylim', [-50,50]); title('PSD of AR(4) Process - Biased');

%% Different figure
figure()
subplot(1, 3, 1); plot(normf, 10*log10(P1_unb)); 
hold on;  plot(normf, 10*log10(P1_b))
xlim([0, 1]); % plotting only for positive frequencies
ylim([-40, 30])
xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14); title('PSD of WGN');
legend('Unbiased', 'Biased')

subplot(1, 3, 2); plot(normf, 10*log10(P2_unb)); 
hold on; plot(normf, 10*log10(P2_b)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14); title('PSD of Noisy Sinusoid');
legend('Unbiased', 'Biased')

subplot(1, 3, 3); plot(normf, 10*log10(P3_unb)); 
hold on; plot(normf, 10*log10(P3_b)); 
xlim([0, 1]); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
set(gca, 'FontSize', 14); title('PSD of AR(4) Process');
legend('Unbiased', 'Biased')
%% Functions
function  r_k = acf(x, biased)

N = length(x);
r_k = zeros(1, N); 
for k=0:N-1
    if strcmp(biased, 'biased')
        a = (1/N); % factor for biased estimate
    else
        a = (1/(N-k)); % factor for unbiased estimate
    end
    r_k(k+1) = a * dot(x(k+1:N), conj(x(1:N-k)));
end

end
function [P, f] = correlogram(x, biased, fs) 
    P = abs(fft(xcorr(x, biased))); 
    L = length(x);
    P = P(1:floor(L/2)+1);
    f = (fs/2)*[0:length(P)- 1]/length(P);
end
