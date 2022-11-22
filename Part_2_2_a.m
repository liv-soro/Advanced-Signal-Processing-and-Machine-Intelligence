% Adaptive Step Sizes
clc; clear; close all

%% a) Gradient Adaptive Step Size (GASS)

% Setup
wo = 0.9; % MA coefficient
order_MA = length(wo)+1;
std_eta = sqrt(0.5); % std of noise
nSamp = 1000; % length of input

eta = std_eta*randn(nSamp, 100); % eta(n)
eta_delayed = zeros(nSamp, 100); 
eta_delayed(2:end, :) = eta(1:end-1, :); % eta(n-1)

x = wo*eta_delayed + eta; % input vector

GASS = {'standard', 'ben', 'af', 'mx'}; % tags for different GASS algorithms
leak = 0;
%%
stepsize = [0.1, 0.01]; % step size / learning rate

% Searching for best hyperparameters
% for alpha = [0.8, 0.9] % alpha = [0.1, 0.5, 0.8] --> [0.8, 0.9]
% for rho = [0.005, 0.01] % rho = [0.01, 0.1, 0.5] --> [0.01, 0.02, 0.001] --> [0.001, 0.005, 0.1]

rho = 0.005; % step-size adaptation parameter
alpha = 0.8; 

% 100 realisations
nIter = 100; % number of iterations
error_tot = zeros(length(stepsize), length(GASS), nIter, nSamp);
weights_tot = zeros(length(stepsize), length(GASS), nIter, nSamp+1, order_MA);

for i = 1:length(stepsize)
    st = stepsize(i);
    for g = 1:length(GASS)
        for iter = 1:nIter
            mdl_in = [eta_delayed(: , iter), eta(:, iter)];
            [weights_tot(i, g, iter, :, :), ~, error_tot(i, g, iter, :)]...
                = LMS_gass(x(:, iter), mdl_in, order_MA, st, rho, leak, GASS{g}, alpha);
        end
    end
end

% Weight Errors
weights_avg = squeeze(mean(weights_tot, 3));
w_e_1_std = wo - squeeze(weights_avg(1,1,:,1));
w_e_2_std = wo - squeeze(weights_avg(2,1,:,1));
w_e_1_ben = wo - squeeze(weights_avg(1,2,:,1));
w_e_1_af = wo - squeeze(weights_avg(1,3,:,1));
w_e_1_mx = wo - squeeze(weights_avg(1,4,:,1));

% Prediction Errors
error_avg = squeeze(mean(error_tot.^2, 3));
e_1_std = squeeze(error_avg(1, 1, :));
e_2_std = squeeze(error_avg(2, 1, :));
e_1_ben = squeeze(error_avg(1, 2, :));
e_1_af = squeeze(error_avg(1, 3, :));
e_1_mx = squeeze(error_avg(1, 4, :));

% Plot
figure()
subplot(1, 3, 2); hold on  
plot(1:nSamp+1, w_e_2_std, 'LineWidth', 1.3);
plot(1:nSamp+1, w_e_1_std, 'LineWidth', 1.3); 

plot(1:nSamp+1, w_e_1_ben, 'LineWidth', 1.3); 
plot(1:nSamp+1, w_e_1_af, 'LineWidth', 1.3); 
plot(1:nSamp+1, w_e_1_mx, 'LineWidth', 1.3); 
xlim([0, nSamp])
xlabel('Time Step')
ylabel('Weight Error')
legend('\mu = 0.01', '\mu = 0.1', 'Benveniste', 'Ang & Farhang', 'Matthews & Xie')
title('GASS Weight Error Curves with \mu_0 = 0.11')
set(gca, 'FontSize', 13);
ylim([-0.02, 0.12]); xlim([0,200])

%%
figure()
subplot(1, 2, 2); 
plot(1:nSamp, 10*log10(e_1_std+eps), 'LineWidth', 1.3); hold on
plot(1:nSamp, 10*log10(e_2_std+eps), 'LineWidth', 1.3);
plot(1:nSamp, pow2db(e_1_ben), 'LineWidth', 1.3);
plot(1:nSamp, 10*log10(e_1_af+eps), 'LineWidth', 1.3);
plot(1:nSamp, 10*log10(e_1_mx+eps), 'LineWidth', 1.3);
xlim([0, nSamp])
legend('\mu = 0.01', '\mu = 0.1', 'Benveniste', 'Ang & Farhang', 'Matthews & Xie')
title("alpha = " + alpha + " and rho = " + rho)


% I AM LOOKING FOR THE ERROR OF THE WEIGHTS!!

% How alpha influences the Ang and Farhang:
%   constant alpha smaller but close to 1 (from their paper)
%   like low pass filter
%   more stable adaptation of step-size parameter than Matthews (alpha = 0)

% How rho influences all algorithms:
%   rho is the learning rate of the learning rate
%   i.e. by how much we value the new update wrt to the past one
%   has to be very small because big jumps in the mearning rate will lead
%   to no convergence

% adv/disadv wrt LMS: complexity up, performance better

sse_std2 = mean(e_2_std(900:end))%0.01
sse_std1 = mean(e_1_std(300:end))
sse_ben = mean(e_1_ben(300:end))
sse_af = mean(e_1_af(300:end))
sse_mx = mean(e_1_mx(300:end))




