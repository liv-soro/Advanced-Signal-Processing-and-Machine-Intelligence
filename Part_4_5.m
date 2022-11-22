%% 4.5. Pretraining Weights
addpath('Data');
load('time-series.mat')

M = 4; % filter order
stepsize = 1e-5; % step size /  learning rate

y_in = delay(y,1); 
[weights, prediction, error] = LMS_tanh(y, y_in, M, stepsize, 39, 1, 1);

% MSE
MSE = mean(error.^2); MSE_dB = 10*log10(MSE);
% Prediction Gain
R = 10*log10(var(prediction)/var(error));

figure(); 
subplot(1, 2, 1); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('True Signal', 'AR(4)- LMS')
title('Time Series AR(4) Prediction with LMS')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('True Signal', 'AR(4)- LMS')
title('Zoomed: AR(4) Prediction with LMS')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);







