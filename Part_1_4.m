% 1.4.a), b) Spectrum Autoregressive Process

% a) based on 1.3.a): shortcomings of unbiased estimate of the autocorrelation 
%   function is that it can be highly erratic when computed for large lags
%   (close to N) where fewer samples are available to estimate the PSD
%   see lecture 3, slide 19 (positive semidefinite)
%   (if variant enough , the estimate might lead to AR parameters of an
%   unstable process)
%   see yule walker equation -> invertibility is not guaranteed for
%   unbiased estimate (and AR param a may not be findable with the inverse
%   Y-W eqn)
% biased form guarantees the matrix being positive definite and hence
% invertible

% b) 
N = [1000, 10000]; %number of samples used to estimate the AR spectral estimate 

a = [2.76, -3.81, 2.65, -0.92]; % coefficients of AR process
a = [1, -a]; % correct for 'filter'

figure();
for i = 1:2
    n = N(i);
    x = filter(1, a, randn(n,1));
    x = x(500:end); % remove transient output of the filter
    
    [H, f] = freqz(1, a, []); % frequency response of the filter
    P = abs(H).^2;
    
    subplot(1, 2, i);
    plot(f/pi, 10*log10(P), 'Linewidth', 1.5, 'Color', [0, 0.4470, 0.7410]); hold on;
    leg = {'True PSD'};
    for order = [4, 8, 10]
        [pxx, w] = pyulear(x, order); % power with Yule-Walker
        plot(w/pi, 10*log10(pxx), 'LineWidth', 1.5)
        title("order: " + order)
        leg{end+1} = sprintf('Order %d', order);
    end
    legend(leg); xlabel('Normalised Frequency'); ylabel('Power (dBs)');
    title("Estimated PSD of AR model (" + N(i) + " samples)")
    set(gca, 'FontSize', 14) 
    hold off
end
