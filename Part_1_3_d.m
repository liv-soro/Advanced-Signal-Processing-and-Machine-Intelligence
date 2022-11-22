% 1.3.d)  Correlation Estimation - Resolution of the Periodogram
clear; close all

points = [30, 60, 120]; 
figure(); hold on
for i = 1:length(points)
    n = points(i);
    N = [0:n];
    res = 1/n
    noise = 0.2/sqrt(2)*(randn(size(N))+ 1i*randn(size(N))); 
    x = exp(1i*2*pi*0.3*N) + exp(1i*2*pi*0.32*N)  + noise;
    acf = xcorr(x, 'biased');
    P = abs(fftshift(fft(acf)));
    normf = linspace(-1, 1, length(P));
    plot(normf, P, 'LineWidth', 1.3)
end
xlim([0,1]); xlabel('Normalised Frequency'); ylabel('Power')
title('PSD of Complex Exponential');
legend('30 samples', '60 samples', '120 samples')
set(gca, 'FontSize', 14);