% 1.3.e)  Correlation Estimation - MUSIC method

% plot for different sample lengths (same as in 1.3 d)
% plot multiple to talk about variance and bias (defo we will have bias)

n = [0:30];
N = length(n); % number of samples
noise = 0.2/sqrt(2)*(randn(size(n))+ 1i*randn(size(n))); 
x = exp(1i*2*pi*0.3*n) + exp(1i*2*pi*0.32*n)  + noise; % noisy singal with two exponential components
fs = 1; % normalised sampling frequency

[X, R] = corrmtx(x, 14, 'mod'); % R is the biased correlation matrix estimated
[S, F] = pmusic(R, 2, [ ], fs, 'corr'); % MUSIC algorithm uses R and the subspace dimention to estimate a pseudospectrum

psd  = abs(fftshift(fft(x)));
f = (-N/2)*(fs/N) : fs/N : (N/2 - 1)*(fs/N);


figure();
subplot(1, 2, 1); plot(F, S, 'LineWidth', 1.5); 
set(gca, 'xlim', [0.25, 0.40], 'FontSize', 14);
xlabel('Frequency (Hz)'); ylabel('Pseudospectrum')
title('Pseudospectrum using MUSIC')

subplot(1, 2, 2); plot(f, psd, 'LineWidth', 1.5); 
set(gca, 'xlim', [0.25, 0.40], 'FontSize', 14);
xlabel('Frequency (Hz)'); ylabel('Periodogram');
title('Periodogram for the same number of points')
% same number of points

% corrmtx 
% Input: 
%   vector 'x' of length 'n' 
%   positive integer 'm' (which is the prediction model order)¨
%   method: matrix computation method
% Output: 
%   X : Data matrix, s.t. conjtransp(X)*X is a biased estimate of
%   the autocorrelation matrix. Returned for autocorrelation matrix
%   estimation. The size of X depends on the matrix computation method. If
%   'mod', X is the 2(n-m)by(m+1) Toepliz matrix that generate an autocorrelation estimate for
%   the length-n data vector x, derived using forward and backward prediction
%   error estimates. Matrix can be used to perform autoregressive parameter
%   estimation using the modified covariance method.
%   R : Biased autocorrelation matrix estimate computed as
%   conjtransp(X)*X. (m+1)by(m+1). By decomposing it we can separate the
%   eigenvectors of the signal from those of noise and estimate the PSD
%   using for example the music method. In other words(?) the matrix can be
%   used to perform autoregressive parameter estimation

% pmusic
% Input:
%   input signal x (in our case a correlation matrix)
%   'p' is the signal subspace dimension (in our case 2 because we have 2
%   complex exponentials.
%   nfft ([]) specifies the integer length of the FFT, nfft, used to 
%   estimate the pseudospectrum.
%   'fs' (1) sample rate in Hz
%   'corr': forces the input argument x to be interpreted as a correlation 
%   matrix rather than matrix of signal data. For this syntax, x must be a 
%   square matrix, and all of its eigenvalues must be nonnegative (i.e. it's positive semidefinite).
% Outputs:
%   S: pseudospectrum estimate of the input signal x
%   F: vector of frequencies (here in Hz bc we provided fs) at which the
%   pseudospectrum is computed

%% Many reps
reps = 25;
S_all = zeros(reps, length(S));
figure()
subplot(1, 2, 1)
for r = 1:reps
    noise = 0.2/sqrt(2)*(randn(size(n))+ 1i*randn(size(n))); 
    x = exp(1i*2*pi*0.3*n) + exp(1i*2*pi*0.32*n)  + noise; % noisy singal with two exponential components
    [X, R] = corrmtx(x, 14, 'mod'); % R is the biased correlation matrix estimated
    [S, F] = pmusic(R, 2, [ ], fs, 'corr'); % MUSIC algorithm uses R and the subspace dimention to estimate a pseudospectrum
    plot(F, S, 'Color', [0.8500 0.3250 0.0980]); hold on
    S_all(r,:) = S;
end
subplot(1, 2, 1); plot(F, mean(S_all), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.3)
set(gca, 'xlim', [0.25, 0.4], 'FontSize', 14)
xlabel('Frequency (Hz)'); ylabel('Pseudospectrum')
title('Pseudospectrum using MUSIC')

subplot(1, 2, 2); plot(f, psd, 'LineWidth', 1.5); 
set(gca, 'xlim', [0.25, 0.40], 'FontSize', 14);
xlabel('Frequency (Hz)'); ylabel('Periodogram');
title('Periodogram for the same number of points')
%% Different assumed model orders
n = [0:30];
N = length(n); % number of samples
fs = 1; % normalised sampling frequency
reps = 10;
S_all = zeros(reps, length(S), 3);

figure();
for p = 1:3
    subplot(1, 3, p); hold on
    for r = 1:reps
        noise = 0.3/sqrt(2)*(randn(size(n))+ 1i*randn(size(n))); 
        x = exp(1i*2*pi*0.3*n) + exp(1i*2*pi*0.32*n)  + noise; % noisy singal with two exponential components
        [X, R] = corrmtx(x, 14, 'mod'); % R is the biased correlation matrix estimated
        [S, F] = pmusic(R, p, [ ], fs, 'corr'); % MUSIC algorithm uses R and the subspace dimention to estimate a pseudospectrum
        plot(F, S, 'Color', [0.9290, 0.6940, 0.1250])
        S_all(r,:,p) = S;
    end
end

subplot(1, 3, 1); plot(F, mean(S_all(:,:,1)), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.3)
set(gca, 'xlim', [0.25, 0.4])
subplot(1, 3, 2); plot(F, mean(S_all(:,:,2)), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.3)
set(gca, 'xlim', [0.25, 0.4])
subplot(1, 3, 3); plot(F, mean(S_all(:,:,3)), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.3)
set(gca, 'xlim', [0.25, 0.4])

% figure();
% subplot(1, 2, 1); plot(F, S, 'LineWidth', 1.5); 
% set(gca, 'xlim', [0.25, 0.40], 'FontSize', 14);
% xlabel('Frequency (Hz)'); ylabel('Pseudospectrum')
% title('Pseudospectrum using MUSIC')
% 
% subplot(1, 2, 2); plot(f, psd, 'LineWidth', 1.5); 
% set(gca, 'xlim', [0.25, 0.40], 'FontSize', 14);
% xlabel('Frequency (Hz)'); ylabel('Periodogram');
% title('Periodogram for the same number of points')
% % same number of points