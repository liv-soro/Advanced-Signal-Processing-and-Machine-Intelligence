% The Least Mean Square Algortihm
clc; clear; close all
%% a) Find the correlation matrix and bounds on the stepsize
% Done on paper
R = [0.926, 0.4629; 0.4629, 0.926];
e = eig(R);

%% b) Implement an LMS adaptive predictor

% Setup
a = [0.1, 0.8]; % coefficients of AR process
order = length(a);
a = [1, -a]; % correct for 'filter'

nSamp = 1000; % length of input

std_eta = sqrt(0.25);
x = filter(1, a, std_eta*randn(nSamp+500, 100)); % 100 realisation of column vector input (AR(2) process)  
x = x(501:end, :); % remove transient effects of the filter

stepsize = [0.05, 0.01]; % step size/ learning rate

% 1 realisation
[weights1, prediction1, error1] = LMS(x(:, 1), order, stepsize(1));
[weights2, prediction2, error2] = LMS(x(:, 1), order, stepsize(2));

sum1 = sum(abs(error1), 'all'); %236.0
sum2 = sum(abs(error2), 'all'); %209.4

% Plot
figure(); subplot(1, 2, 1); hold on;
plot(1:nSamp, 10*log10(error1.^2));
plot(1:nSamp, 10*log10(error2.^2));
leg = legend('\mu = 0.05', '\mu = 0.01');
xlabel('Time Step'); ylabel('Squared Error (dBs)');
title('Squared Prediction Error (1 realisation)')
set(gca, 'FontSize', 14) 
% ylim([-80, 0])
hold off

% 100 realisations
nIter = 100; % number of iterations
error_tot = zeros(2, nIter, nSamp);
for i = 1:2
    st = stepsize(i);
    for iter = 1:nIter
        [weights, prediction, error_tot(i, iter, :)] = LMS(x(:, iter), order, st);
    end
end

error_tot = squeeze(mean(10*log10(error_tot.^2), 2)); % average over realisation

% Plot
subplot(1, 2, 2); hold on
plot(1:nSamp, error_tot(1, :));
plot(1:nSamp, error_tot(2, :));
% ylim([-25, -10])
leg = legend('\mu = 0.05', '\mu = 0.01');
xlabel('Time Step'); ylabel('Mean Squared Error (dBs)');
title('Learning Curve (100 realisations)')
set(gca, 'FontSize', 14) 
hold off
% title(leg, 'Step Size')
%% c) Misadjustment
clear; clc;
R = [0.926, 0.4629; 0.4629, 0.926];
T = trace(R);

% Setup
a = [0.1, 0.8]; % coefficients of AR process
order = length(a);
a = [1, -a]; % correct for 'filter'

nSamp = 100000; % length of input
stepsize = [0.01, 0.05]; % step size/ learning rate
std_eta = sqrt(0.25);

x = filter(1, a, std_eta*randn(nSamp+500, 100)); % 100 realisation of column vector input (AR(2) process)  
x = x(501:end, :); % remove transient effects of the filter

% 100 realisations
nIter = 100; % number of iterations
error_tot = zeros(2, nIter, nSamp);
for i = 1:2
    st = stepsize(i);
    for iter = 1:nIter
        [~, ~, error_tot(i, iter, :)] = LMS(x(:, iter), order, st);
    end
end
MSE = squeeze(mean(error_tot(:, :, 90001:end).^2, 2)); % average over realisation
MSE = mean(MSE, 2); 

% Calculate Estimated Misadjustment (MSE = std_eta + EMSE; Mis = EMSE/std_eta)
Mis1_est = (MSE(1) - std_eta.^2) / std_eta.^2 ; % 0.0029 % 0.0074
Mis2_est = (MSE(2) - std_eta.^2) / std_eta.^2 ; % 0.0127 % 0.0516

% Calculate Theoretical Misadjustment
Mis1_th = (stepsize(1)/2)*T; % 0.0093 % 0.0093
Mis2_th = (stepsize(2)/2)*T; % 0.0463 % 0.0463

percentage_diff1 = abs(Mis1_th - Mis1_est)/Mis1_th; % 0.2047
percentage_diff2 = abs(Mis2_th - Mis2_est)/Mis2_th; % 0.1154

% Mis1_th is smaller than Mis2_th because smaller step size
% Theres a bigger error for bigger step size because the theoretical value
% is less well approximated by (mu/2)*T

%% d) Steady State Filter Coefficients
clear; clc;
% Setup
a = [0.1, 0.8]; % coefficients of AR process
order = length(a);
a = [1, -a]; % correct for 'filter'

nSamp = 5000; % length of input
stepsize = [0.05, 0.01]; % step size/ learning rate
std_eta = sqrt(0.25);

x = filter(1, a, std_eta*randn(nSamp+500, 100)); % 100 realisation of column vector input (AR(2) process)  
x = x(501:end, :); % remove transient effects of the filter

% 100 realisations
nIter = 100; % number of iterations
weights_tot = zeros(2, nIter, nSamp+1, order);
for i = 1:2
    st = stepsize(i);
    for iter = 1:nIter
        [weights_tot(i, iter, :, :), ~, ~] = LMS(x(:, iter), order, st);
    end
end

mean_weights_1 = squeeze(mean(weights_tot(1, :, :, :), 2)); % 0.0948, 0.7919
mean_weights_2 = squeeze(mean(weights_tot(2, :, :, :), 2)); % 0.0819, 0.7685

%% Plot
figure(); subplot(1, 2, 1); hold on
plot(1:nSamp+1, mean_weights_1(:, 1), 'LineWidth', 1.3, 'Color', [0, 0.4470, 0.7410])
plot(1:nSamp+1, mean_weights_1(:, 2), 'LineWidth', 1.3, 'Color', [0.8500, 0.3250, 0.0980])
plot(1:nSamp+1, 0.1*ones(1, nSamp+1), '--','LineWidth', 1.3, 'Color', [0, 0.4470, 0.7410])
plot(1:nSamp+1, 0.8*ones(1, nSamp+1), '--', 'LineWidth', 1.3, 'Color', [0.8500, 0.3250, 0.0980])

legend('â_1', 'â_2', 'a_1', 'a_2'); 
xlabel('Time Step'); ylabel('Parameter Estimate');
title({'Filter Coefficients Estimates', 'with LMS (\mu = 0.05)'})
set(gca, 'FontSize', 14) 
xlim([0, 5000]); ylim([0, 0.9])
hold off

subplot(1, 2, 2); hold on
plot(1:nSamp+1, mean_weights_2(:, 1), 'LineWidth', 1.3, 'Color', [0, 0.4470, 0.7410])
plot(1:nSamp+1, mean_weights_2(:, 2), 'LineWidth', 1.3, 'Color', [0.8500, 0.3250, 0.0980])
plot(1:nSamp+1, 0.1*ones(1, nSamp+1), '--','LineWidth', 1.3, 'Color', [0, 0.4470, 0.7410])
plot(1:nSamp+1, 0.8*ones(1, nSamp+1), '--', 'LineWidth', 1.3, 'Color', [0.8500, 0.3250, 0.0980])

legend('â_1', 'â_2', 'a_1', 'a_2'); 
xlabel('Time Step'); ylabel('Parameter Estimate');
title({'Filter Coefficients Estimates', 'with LMS (\mu = 0.01)'})
set(gca, 'FontSize', 14) 
xlim([0, 5000]); ylim([0, 0.9])
hold off

% for smaller step size
% slower convergence
% less error

% 0.05
a1_mu1 = mean(mean_weights_1(end-100:end, 1)); % 0.0716
a2_mu1 = mean(mean_weights_1(end-100:end, 2)); % 0.7207
% 0.01
a1_mu2 = mean(mean_weights_2(end-100:end, 1)); % 0.0878
a2_mu2 = mean(mean_weights_2(end-100:end, 2)); % 0.7766

%% f) Leaky LMS

% Setup
a = [0.1, 0.8]; % coefficients of AR process
order = length(a);
a = [1, -a]; % correct for 'filter'

nSamp = 1500; % length of input

x = filter(1, a, 0.5*randn(nSamp+500, 100)); % 100 realisation of column vector input (AR(2) process)  
x = x(501:end, :); % remove transient effects of the filter

stepsize = [0.05, 0.01]; % step size/ learning rate
leakage = [0.05, 0.1, 0.5]; % leaky coefficients

% % 1 realisation
% [weights1, prediction1, error1] = LMS_leaky(x(:, 1), order, stepsize(1), leakage(1));
% [weights2, prediction2, error2] = LMS_leaky(x(:, 1), order, stepsize(2), leakage(1));

% 100 realisations
nIter = 100; % number of iterations
weights_tot = zeros(2, 3, nIter, nSamp+1, order);
for i = 1:2
    st = stepsize(i);
    for j = 1:3
        leak = leakage(j);
        for iter = 1:nIter
            [weights_tot(i, j, iter, :, :), ~, ~] = LMS_leaky(x(:, iter), order, st, leak);
        end
    end
end

w_avg = cell(2, 3);
for i=1:2
    for j=1:3
        w_avg{i,j} = squeeze(mean(weights_tot(i, j, :, :, :), 3)); 
    end
end

% Subplot
figure();
sp = 0;
for j=1:3
    for i=1:2
        sp = sp+1;
        subplot(3, 2, sp); hold on
        plot(1:nSamp+1, w_avg{i,j}(:,1), 'LineWidth', 1.3, 'Color', [0, 0.4470, 0.7410])
        plot(1:nSamp+1, w_avg{i,j}(:,2), 'LineWidth', 1.3, 'Color', [0.8500, 0.3250, 0.0980])
        plot(1:nSamp+1, 0.1*ones(1, nSamp+1), '--', 'LineWidth', 1.3, 'Color', [0, 0.4470, 0.7410])
        plot(1:nSamp+1, 0.8*ones(1, nSamp+1), '--', 'LineWidth', 1.3, 'Color', [0.8500, 0.3250, 0.0980])
        legend('â_1', 'â_2', 'a_1', 'a_2')
        ylim([0,1]); xlim([0, 1500])
        title("Estimates with \mu = " + stepsize(i) + " and \gamma = " + leakage(j))
        ylabel('Estimate'); xlabel('Time Step')
        set(gca, 'FontSize', 12)
        hold off
    end
end













