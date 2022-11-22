%% 4. From LMS to Deep Learning
%% 4.1. Time-Series Prediction with Standard LMS
addpath('Data');
load('time-series.mat')

mean_y = mean(y);
y = detrend(y - mean_y); % remove mean and detrend

M = 4; % filter order
stepsize = 1e-5; % step size /  learning rate

% LMS
[weights, prediction, error] = LMS(y, M , stepsize);

% MSE
MSE = mean(error.^2); MSE_dB = 10*log10(MSE);

% Prediction Gain
R = 10*log10(var(prediction)/var(error));

% Zero-mean version of y against its one-step ahead prediction
figure(); 
subplot(1, 2, 1); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('True Signal', 'AR(4)- LMS')
title('AR(4) Prediction on Time Series with LMS')
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


% This algorithm performs well for a linear model - unable to capture
% nonlinearities within the signal. Fails to perfectly adjusts to the
% original signal.

%% 4.2. Dynamical Perceptron
y_in = delay(y, 1);

% LMS with activation function
[weights, prediction, error] = LMS_tanh(y, y_in, M , stepsize, 1, 0, 0);

% MSE
MSE = mean(error.^2); MSE_dB = 10*log10(MSE);

% Prediction Gain
R = 10*log10(var(prediction)/var(error));

% Zero-mean version of y against its one-step ahead prediction
figure(); 
subplot(1, 2, 1); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('True Signal', 'LMS - tanh')
title('Prediction with Dynamical Perceptron (tanh)')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('True Signal', 'LMS - tanh')
title('Prediction with Dynamical Perceptron (tanh)')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);


%% 4.3. Dynamical Perceptron with Scaled Activation Function
scales = [1, 10:10:100]; % scaling of activation function
leak = 0;

y_in = delay(y, M);
MSE = zeros(1, length(scales));
R = zeros(1, length(scales));

for i = 1:length(scales)
    % LMS with activation function
    [weights, prediction, error] = LMS_tanh(y_in, y, M , stepsize, leak, scales(i));

    % MSE
    MSE(i) = mean(error.^2); MSE_dB = 10*log10(MSE);

    % Prediction Gain
    R(i) = 10*log10(var(prediction)/var(error));
end

figure(); hold on
yyaxis left;
plot(scales, MSE_dB, 'Linewidth', 1.4)
ylabel('MSE (dB)')
yyaxis right;
plot(scales, R,  'Linewidth', 1.4)
ylabel('Prediction Gain, R')
xlabel('Scale of the Activation Function')

title('MSE and R vs Scaling of the Activation Function')
set(gca, 'FontSize', 13)
%% Zero-mean version of y against its one-step ahead prediction
[weights, prediction, error] = LMS_tanh(y_in, y, M , stepsize, leak, 39);
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


