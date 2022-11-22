%% 4. From LMS to Deep Learning
clear; close all

%%
addpath('Data');
load('time-series.mat')
MSE = zeros(1, 5); % Mean squared error
R = zeros(1, 5); % Prediction Gain
intv = 1:250;
%% 4.1) Time-Series Prediction with Standard LMS

% Remove mean and detrend
mean_y = mean(y);
y_zeroavg = detrend(y - mean_y); 

M = 4; % filter order
stepsize = 1e-5; % step size /  learning rate

% LMS
[~, prediction, error] = LMS(y_zeroavg, M , stepsize);
MSE(1) = mean(error.^2); 
R(1) = 10*log10(var(prediction)/var(error));

figure(); subplot(1, 2, 1); hold on
plot(1:length(y), y_zeroavg, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('AR(4) Prediction on zero-avg Time Series')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y_zeroavg, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Segment')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

% This algorithm performs well for a linear model - unable to capture
% nonlinearities within the signal. Fails to perfectly adjusts to the
% original signal.

%% 4.2. Dynamical Perceptron
y_in = delay(y_zeroavg, 1);

% LMS with activation function
[~, prediction, error] = LMS_tanh(y_zeroavg, y_in, M , stepsize, 1, 0, 0);

MSE(2) = mean(error.^2); 
R(2) = 10*log10(var(prediction)/var(error));

figure(); subplot(1, 2, 1); hold on
plot(1:length(y), y_zeroavg, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Prediction with Activation Function (tanh)')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y_zeroavg, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Segment')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

%% 4.3. Dynamical Perceptron with Scaled Activation Function
MSE_3 = zeros(1, length(intv));
R_3 = zeros(1, length(intv));
count = 0;
for i = intv
   count = count+1;
   [~, prediction, error] = LMS_tanh(y_zeroavg, y_in, M , stepsize, i, 0, 0); 
   MSE_3(count) = mean(error.^2); 
   R_3(count) = 10*log10(var(prediction)/var(error));
end
figure(); 
yyaxis left; plot(intv, 10*log10(MSE_3), 'LineWidth', 1.4); ylabel('MSE (dB)')
hold on;
yyaxis right; plot(intv, R_3, 'LineWidth', 1.4); ylabel('Prediction Gain (dB)')
xlabel('Scale of the Activation Function')
title('MSE and Prediction Gain vs Scaling')
set(gca, 'FontSize', 13)

% LMS with scaled activation function on zero-mean time series
[~, prediction, error] = LMS_tanh(y_zeroavg, y_in, M , stepsize, 82, 0, 0);

MSE(3) = mean(error.^2); 
R(3) = 10*log10(var(prediction)/var(error));

figure(); subplot(1, 2, 1); hold on
plot(1:length(y), y_zeroavg, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Prediction with Scaled Activation Function (tanh)')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y_zeroavg, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Segment')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

%% 4.4. Dynamical Perceptron with Bias
y_in = delay(y, 1);

MSE_4 = zeros(1, length(intv));
R_4 = zeros(1, length(intv));
count = 0;
for i = intv
   count = count+1;
   [~, prediction, error] = LMS_tanh(y, y_in, M , stepsize, i, 1, 0); 
   MSE_4(count) = mean(error.^2); 
   R_4(count) = 10*log10(var(prediction)/var(error));
end
figure(); 
yyaxis left; plot(intv, 10*log10(MSE_4), 'LineWidth', 1.4); ylabel('MSE (dB)')
hold on;
yyaxis right; plot(intv, R_4, 'LineWidth', 1.4); ylabel('Prediction Gain (dB)')
xlabel('Scale of the Activation Function')
title('MSE and Prediction Gain vs Scaling')
set(gca, 'FontSize', 13)

% LMS with bias on original time series
[~, prediction, error] = LMS_tanh(y, y_in, M, stepsize, 53, 1, 0);

MSE(4) = mean(error.^2); 
R(4) = 10*log10(var(prediction)/var(error));

figure(); subplot(1, 2, 1); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Prediction with bias')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Segment')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

%% 4.5. Pretraining Weights

MSE_5 = zeros(1, length(intv));
R_5 = zeros(1, length(intv));
count = 0;
for i = intv
   count = count+1;
   [~, prediction, error] = LMS_tanh(y, y_in, M , stepsize, i, 1, 1); 
   MSE_5(count) = mean(error.^2); 
   R_5(count) = 10*log10(var(prediction)/var(error));
end
figure(); 
yyaxis left; plot(intv, 10*log10(MSE_5), 'LineWidth', 1.4); ylabel('MSE (dB)')
hold on;
yyaxis right; plot(intv, R_5, 'LineWidth', 1.4); ylabel('Prediction Gain (dB)')
xlabel('Scale of the Activation Function')
title('MSE and Prediction Gain vs Scaling')
set(gca, 'FontSize', 13)

% Using pre trained weights to do LMS with bias on original time series
[~, prediction, error] = LMS_tanh(y, y_in, M, stepsize, 194, 1, 1);

MSE(5) = mean(error.^2); 
R(5) = 10*log10(var(prediction)/var(error));

figure(); subplot(1, 2, 1); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Prediction with Pretraining of Weights')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

subplot(1, 2, 2); hold on
plot(1:length(y), y, 'LineWidth', 1.3); 
plot(1:length(prediction), prediction, 'LineWidth', 1.3)
xlim([800,1000])
legend('$y$', '$\hat{y}$', 'Interpreter', 'Latex')
title('Segment')
xlabel('Time Step (n)'); ylabel('Amplitude')
set(gca, 'FontSize', 13);

%% MSE and Prediction Gain vs Scaling
figure(); 
subplot(1, 3, 1); hold on
yyaxis left; plot(intv, 10*log10(MSE_3), 'LineWidth', 1.4); ylabel('MSE (dB)')
hold on;
yyaxis right; plot(intv, R_3, 'LineWidth', 1.4); ylabel('Prediction Gain (dB)')
xlabel('Scale')
title('4.3')
set(gca, 'FontSize', 13)
xlim([1,250])

subplot(1, 3, 2); hold on
yyaxis left; plot(intv, 10*log10(MSE_4), 'LineWidth', 1.4); ylabel('MSE (dB)')
hold on;
yyaxis right; plot(intv, R_4, 'LineWidth', 1.4); ylabel('Prediction Gain (dB)')
xlabel('Scale')
title('4.4')
set(gca, 'FontSize', 13)
xlim([1,250])

subplot(1, 3, 3); hold on
yyaxis left; plot(intv, 10*log10(MSE_5), 'LineWidth', 1.4); ylabel('MSE (dB)')
hold on;
yyaxis right; plot(intv, R_5, 'LineWidth', 1.4); ylabel('Prediction Gain (dB)')
xlabel('Scale')
title('4.5')
set(gca, 'FontSize', 13)
xlim([1,250])

sgtitle('MSE and Prediction Gain vs Scale of the Activation Function', 'FontSize', 16)