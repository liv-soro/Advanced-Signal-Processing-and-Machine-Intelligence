% 3. Widely Liner Filtering and Adaptive Spectrum Estimation
clear; close all; clc

%% Complex LMS and Widely Linear Modelling
% a) Learning curves of CLMS and ACLMS for a WLMA(1) process

% Generate first order widely linear moving average process (WLMA)
nIter = 100; 
x = wgn(nIter, 1000, 0, 'complex'); % circular WGN
b1 = 1.5 + 1i;  b2 = 2.5 - 0.5i; % coefficients
y = x + b1 * [zeros(nIter, 1), x(:, 1:end-1)] + b2 * [zeros(nIter, 1), conj(x(:, 1:end-1))]; % WLMA
stepsize = 0.1;
M = 2; 

e_CLMS = zeros(1000, nIter);
e_ACLMS = zeros(1000, nIter);

% Learning curves of CLMS and ACLMS on a WLMA process of order 1 
for iter = 1:nIter
    y_delay = [0, y(iter, 1:end-1)];
    x_in = delay(x(iter, :)', M);
    [~, ~, e_CLMS(:, iter)] = CLMS(x_in', y_delay, stepsize, M);
    [~, ~, ~, e_ACLMS(:, iter)] = ACLMS(x_in', y_delay, stepsize, M);
end

e_avg_CLMS = mean(abs(e_CLMS).^2, 2);
e_avg_ACLMS = mean(abs(e_ACLMS).^2, 2);

%%
figure(); subplot(2, 1, 1); hold on
plot(1:1000, 10*log10(e_avg_CLMS), 'LineWidth', 1.3)
plot(1:1000, 10*log10(e_avg_ACLMS), 'LineWidth', 1.3)
legend('CLMS', 'ACLMS')
title({'Learning Curves for a WLMA(1) of CLMS and ACLMS'}, 'FontSize', 13)
xlabel('Time Step (n)', 'FontSize', 13)
ylabel('Mean Square Error (dBs)', 'FontSize', 13)
% as expected, CLMS doesnt learn. This is bc it is a generic extension of its 
% real counterpart and therefore it works only for circular random
% variables (WLMA isnt, as seen below)

% Circularity of the complex vectors:
[~, circ_coeffWGN] = circularity(x);
[~, circ_coeffWLMA] = circularity(y);

subplot(2, 2, 3); hold on;
scatter(real(x(:)), imag(x(:)))
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title('White Noise |\rho| = 0.12', 'FontSize', 13)

subplot(2, 2, 4); hold on;
scatter(real(y(:)), imag(y(:)), [], [0.8500, 0.3250, 0.0980])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title('WLMA(1) |\rho| = 0.86', 'FontSize', 13)

