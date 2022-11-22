%% c) Generalized Normalised Gradient Descent

% Setup
wo = 0.9; % MA coefficient
order_MA = length(wo)+1;
std_eta = sqrt(0.5); % std of noise
nSamp = 1000; % length of input

eta = std_eta*randn(nSamp, 100); % eta(n)
eta_delayed = zeros(nSamp, 100); 
eta_delayed(2:end, :) = eta(1:end-1, :); % eta(n-1)

x = wo*eta_delayed + eta; % input vector

GASS = 'ben'; % tags for different GASS algorithms
leak = 0;
stepsize = 0.1 %[0.01, 0.1, 1]; % step size / learning rate
rho = 0.005; 
alpha = NaN; 

% 100 realisations
nIter = 100; % number of iterations
weights_ben = zeros(length(stepsize), nIter, nSamp+1, order_MA);
error_ben = zeros(length(stepsize), nIter, nSamp);
weights_GNGD = zeros(length(stepsize), nIter, nSamp+1, order_MA);

ben_w_avg = zeros(length(stepsize), nSamp+1, order_MA);
GNGD_w_avg = zeros(length(stepsize), nSamp+1, order_MA);


for i = 1:length(stepsize)
    st = stepsize(i);
    for iter = 1:nIter
        mdl_in = [eta_delayed(: , iter), eta(:, iter)];
        [weights_ben(i, iter, :, :), ~, error_ben(i, iter, :)] ...
            = LMS_gass(x(:, iter), mdl_in, order_MA, st, rho, leak, GASS, alpha);
        [weights_GNGD(i, iter, :, :), ~, ~] ...
            = GNGD(x(:, iter), mdl_in, order_MA, st, rho);
    end
    ben_w_avg(i, :, :) = squeeze(mean(weights_ben(i, :, :, :), 2));
    GNGD_w_avg(i, :, :) = squeeze(mean(weights_GNGD(i, :, :, :), 2));
    figure();
    plot(1:nSamp+1, wo - ben_w_avg(i, :, 1))
    plot(1:nSamp+1, wo - GNGD_w_avg(i, :, 1))
end

% Weight Errors
%figure();
%ben_avg_1 = squeeze(mean(weights_ben(1,:,:,:), 2));

%% c) With 'optimal parameters'

% Setup
wo = 0.9; % MA coefficient
order_MA = length(wo)+1;
std_eta = sqrt(0.5); % std of noise
nSamp = 1000; % length of input

eta = std_eta*randn(nSamp, 100); % eta(n)
eta_delayed = zeros(nSamp, 100); 
eta_delayed(2:end, :) = eta(1:end-1, :); % eta(n-1)

x = wo*eta_delayed + eta; % input vector

GASS = 'ben'; % tags for different GASS algorithms
leak = 0;
mu_GNGD = 0.1; % step size / learning rate
rho_GNGD = 0.005; 
mu_ben = 0.1;
rho_ben = 0.005;
alpha = NaN; 

% 100 realisations
nIter = 100; % number of iterations
weights_ben = zeros(nIter, nSamp+1, order_MA);
error_ben = ones(nIter, nSamp);
weights_GNGD = zeros(nIter, nSamp+1, order_MA);
error_GNGD = ones(nIter, nSamp);

ben_w_avg = zeros(nSamp+1, order_MA);
GNGD_w_avg = zeros(nSamp+1, order_MA);

for iter = 1:nIter
    mdl_in = [eta_delayed(: , iter), eta(:, iter)];
    [weights_ben(iter, :, :), ~, error_ben(iter, :)] ...
        = LMS_gass(x(:, iter), mdl_in, order_MA, mu_ben, rho_ben, leak, GASS, alpha);
    [weights_GNGD(iter, :, :), ~, error_GNGD(iter, :)] ...
        = GNGD(x(:, iter), mdl_in, order_MA, mu_GNGD, rho_GNGD);
end

ben_w_avg(:, :) = squeeze(mean(weights_ben(:, :, :), 1));
GNGD_w_avg(:, :) = squeeze(mean(weights_GNGD(:, :, :), 1));

ben_e_avg = mean(mag2db((error_ben + eps).^2), 1);
GNGD_e_avg = mean(mag2db((error_GNGD + eps).^2), 1);


%% Plot
figure();
subplot(1, 2, 1)
plot(1:nSamp+1, wo - ben_w_avg(:, 1), 'LineWidth', 1.3); hold on
plot(1:nSamp+1, wo - GNGD_w_avg(:, 1), 'LineWidth', 1.3)
legend('Benveniste', 'GNGD');
xlabel('Time Step'); ylabel('Weight Error')
xlim([0, 80]); ylim([-0.01, 0.91])
title('Weight Error Curves (\mu_0 = 0.1, \rho = 0.005)')
set(gca, 'FontSize', 13);

subplot(1, 2, 2);
plot(1:nSamp, ben_e_avg, 'LineWidth', 1.3); hold on
plot(1:nSamp, GNGD_e_avg, 'LineWidth', 1.3); hold on
legend('Benveniste', 'GNGD');
xlabel('Time Step'); ylabel('Prediciton Error')
%xlim([0, 80]); ylim([-0.01, 0.91])
title('Squared Prediction Error Curves (\mu_0 = 0.1, \rho = 0.005)')
set(gca, 'FontSize', 13);

% Weight Errors
%figure();
%ben_avg_1 = squeeze(mean(weights_ben(1,:,:,:), 2));
