load('PCAPCR.mat')

%% 1.6 a)
% Singular value decomposition of X and Xnoise
[Ux, Sx, Vx] = svd(X);
[Uxn, Sxn, Vxn] = svd(Xnoise);

% Plotting Singular Values
figure(); subplot(1, 2, 1); hold on;
stem(diag(Sx), 'LineWidth', 1.4); stem(diag(Sxn), '--x', 'LineWidth', 1.4); xlim([0,10]);
set(gca, 'Fontsize', 13); title('Singular Values of Input and Noisy Input Signals')
xlabel('Subspace Dimension Index'); ylabel('Singular Value');
legend('Input', 'Noisy Input')

% rank of input data = 3 (observation)
r = rank(X);

% Square error
subplot(1, 2, 2); 
stem((diag(Sx)-diag(Sxn)).^2, 'LineWidth', 1.4); xlim([0, 10])
set(gca, 'Fontsize', 13); title({'Square Error between Singular Values', 'of Input and Noisy Input Singals'})
xlabel('Subspace Dimension Index'); ylabel('Singular Value Square Error');


%% 1.6 b) Denoising
% try other way of doing it 
[Uxn_dash, Sxn_dash, Vxn_dash] = svds(Xnoise,3); % 3 largest singular values
Xn_dash = Uxn_dash*(Sxn_dash*Vxn_dash'); % denoised input

% Comparing error between variables of input X with Xnoise and Xdenoised
% (X_dash)
figure(); hold on;
stem(mean((X-Xnoise).^2,1), 'LineWidth', 1.4); 
stem(mean((X-Xn_dash).^2,1), '--x', 'LineWidth', 1.4);
set(gca, 'Fontsize', 13); title('Error between Signals and Input')
xlabel('Variables'); ylabel('Mean Square Error');
legend('Noisy Input', 'Denoised Input'); ylim([0, 0.33])

e_noise = mean((X-Xnoise).^2, 1);
e_denoise = mean((X-Xn_dash).^2,1);
%% 1.6 c)
% Coefficients for regression matrix using Ordinary Least Squares
B_ols = (transpose(Xnoise) * Xnoise) \ transpose(X) * Y;

% Coefficients for regression matrix using PCA and retaining the r largest
% principal components
B_pcr = Vxn(:, 1:r) * inv(Sxn(1:r, 1:r)) * transpose(Uxn(:, 1:r)) * Y;

% Training output by OLS regression model & PCR regression model
Y_ols = Xnoise * B_ols;
Y_pcr = Xn_dash * B_pcr;

% Test output by OLS regression model
Ytest_ols = Xtest * B_ols;

% Test output by PCR regression model
[Uxt_dash, Sxt_dash, Vxt_dash] = svds(Xtest,3); % 3 largest singular values
Xt_dash = Uxt_dash*(Sxt_dash*Vxt_dash');
Ytest_pcr = Xt_dash * B_pcr;


% Estimation Errors
error_est_ols = mean((Y - Y_ols).^2, 1);
error_est_pcr = mean((Y - Y_pcr).^2, 1);

figure(); 
subplot(1, 2, 1); hold on;
stem(error_est_ols, 'LineWidth', 1.4); stem(error_est_pcr, '--x', 'LineWidth', 1.4);
set(gca, 'Fontsize', 13); title('Estimation Error')
xlabel('Variables'); ylabel('Mean Square Error');
legend('OLS', 'PCR'); 

% Test Errors
error_test_ols = mean((Ytest - Ytest_ols).^2, 1);
error_test_pcr = mean((Ytest - Ytest_pcr).^2, 1);

subplot(1, 2, 2); hold on;
stem(error_test_ols, 'LineWidth', 1.4); 
stem(error_test_pcr, '--x', 'LineWidth', 1.4);
set(gca, 'Fontsize', 13); title('Test Error')
xlabel('Variables'); ylabel('Mean Square Error');
legend('OLS', 'PCR'); 

%% 1.6 d)
reps = 100; 
error_ols = cell(reps, 1);
error_pcr = cell(reps, 1);

for i = 1:reps
    % OLS
    [Y_ols_hat, Y_ols] = regval(B_ols);
    error_ols{i, 1} = mean((Y_ols - Y_ols_hat).^2, 1);
    % PCR
    [Y_pcr_hat, Y_pcr] = regval(B_pcr);
    error_pcr{i, 1} = mean((Y_pcr - Y_pcr_hat).^2, 1);
end

avg_error_ols = mean(cell2mat(error_ols), 1);
avg_error_pcr = mean(cell2mat(error_pcr), 1);

figure(); hold on
stem(avg_error_ols, 'LineWidth', 1.4);
stem(avg_error_pcr, '--x', 'LineWidth', 1.4)
set(gca, 'Fontsize', 13); title('Error between Estimates and Data using regval')
xlabel('Variables'); ylabel('Mean Square Error');
legend('OLS', 'PCR'); 








