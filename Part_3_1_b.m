% 3. Widely Liner Filtering and Adaptive Spectrum Estimation
clear; close all; clc

%% 3.1. Complex LMS and Widely Linear Modelling
% b) CLMS and ACLMS on bivariate Wind Data + exploring effect of filter
% order
high = load('Data/high-wind.mat');
medium = load('Data/medium-wind.mat');
low = load('Data/low-wind.mat');

v_high = high.v_east + high.v_north * 1i;
v_medium = medium.v_east + medium.v_north * 1i;
v_low = low.v_east + low.v_north * 1i;

% Scatter Plots (Circularity Plots)
figure(); subplot(1, 3, 1); 
[~, circ_coeffHI] = circularity(v_high);
scatter(real(v_high(:)), imag(v_high(:)), [], [0.6350, 0.0780, 0.1840])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title("High Wind, |\rho| = " + circ_coeffHI, 'FontSize', 13)
axis equal
mi_re = min(real(v_high(:))); ma_re = max(real(v_high(:)));
mi_im = min(imag(v_high(:))); ma_im = max(imag(v_high(:)));
mi = min([mi_re, mi_im]); ma = max([ma_re, ma_im]);
xlim([mi - abs(mi/3), ma + abs(ma/3)])
ylim([mi - abs(mi/3), ma + abs(ma/3)])

subplot(1, 3, 2); 
[~, circ_coeffME] = circularity(v_medium);
scatter(real(v_medium(:)), imag(v_medium(:)), [], [0.8500, 0.3250, 0.0980])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title("Medium Wind, |\rho| = " + circ_coeffME, 'FontSize', 13)
axis equal
% mi_re = min(real(v_medium(:))); ma_re = max(real(v_medium(:)));
% mi_im = min(imag(v_medium(:))); ma_im = max(imag(v_medium(:)));
xlim([mi - abs(mi/3), ma + abs(ma/3)])
ylim([mi - abs(mi/3), ma + abs(ma/3)])

subplot(1, 3, 3);
[~, circ_coeffLO] = circularity(v_low);
scatter(real(v_low(:)), imag(v_low(:)), [], [0.9290, 0.6940, 0.1250])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title("Low Wind, |\rho| = " + circ_coeffLO, 'FontSize', 13)
axis equal
% mi_re = min(real(v_low(:))); ma_re = max(real(v_low(:)));
% mi_im = min(imag(v_low(:))); ma_im = max(imag(v_low(:)));
xlim([mi - abs(mi/3), ma + abs(ma/3)])
ylim([mi - abs(mi/3), ma + abs(ma/3)])

%% OTHER
% Scatter Plots (Circularity Plots)
figure(); subplot(1, 3, 1); 
[~, circ_coeffHI] = circularity(v_high);
scatter(real(v_high(:)), imag(v_high(:)), [], [0.6350, 0.0780, 0.1840])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title("High Wind, |\rho| = " + circ_coeffHI, 'FontSize', 13)
axis equal
mi_re = min(real(v_high(:))); ma_re = max(real(v_high(:)));
mi_im = min(imag(v_high(:))); ma_im = max(imag(v_high(:)));
mi = min([mi_re, mi_im]); ma = max([ma_re, ma_im]);
xlim([mi - abs(mi/3), ma + abs(ma/3)])
ylim([mi - abs(mi/3), ma + abs(ma/3)])

subplot(1, 3, 2); 
[~, circ_coeffME] = circularity(v_medium);
scatter(real(v_medium(:)), imag(v_medium(:)), [], [0.8500, 0.3250, 0.0980])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title("Medium Wind, |\rho| = " + circ_coeffME, 'FontSize', 13)
axis equal
mi_re = min(real(v_medium(:))); ma_re = max(real(v_medium(:)));
mi_im = min(imag(v_medium(:))); ma_im = max(imag(v_medium(:)));
mi = min([mi_re, mi_im]); ma = max([ma_re, ma_im]);
xlim([mi - abs(mi/3), ma + abs(ma/3)])
ylim([mi - abs(mi/3), ma + abs(ma/3)])

subplot(1, 3, 3);
[~, circ_coeffLO] = circularity(v_low);
scatter(real(v_low(:)), imag(v_low(:)), [], [0.9290, 0.6940, 0.1250])
xlabel('Real Part', 'FontSize', 13);
ylabel('Imaginary Part', 'FontSize', 13);
title("Low Wind, |\rho| = " + circ_coeffLO, 'FontSize', 13)
axis equal
mi_re = min(real(v_low(:))); ma_re = max(real(v_low(:)));
mi_im = min(imag(v_low(:))); ma_im = max(imag(v_low(:)));
mi = min([mi_re, mi_im]); ma = max([ma_re, ma_im]);
xlim([mi - abs(mi/3), ma + abs(ma/3)])
ylim([mi - abs(mi/3), ma + abs(ma/3)])

%% Error Curves against Filter Order
clear, close all; clc
%%
load('Data/high-wind.mat');
wind(1,:) = (v_east + v_north * 1i).';
load('Data/medium-wind.mat');
wind(2,:) = (v_east + v_north * 1i).';
load('Data/low-wind.mat');
wind(3,:) = (v_east + v_north * 1i).';

nWinds = 3; % number of winds
M = 1:24; % filter orders
%stepsize = [0.001, 0.001, 0.001]; % learning rate
stepsize = [0.001, 0.01, 0.1]; % learning rate

e_CLMS = cell(nWinds, length(M));
e_ACLMS = cell(nWinds, length(M));
se_CLMS = zeros(nWinds, length(M));
se_ACLMS = zeros(nWinds, length(M));


for w = 1:nWinds
    st = w;
    for i = 1:length(M)
        [v_in] = delay(wind(w,:)', M(i)); % delay the signal by 1 to filter length
        [~, ~, e_CLMS{w,i}] = CLMS(v_in', wind(w, :), stepsize(st), M(i));
        [~, ~, ~, e_ACLMS{w,i}] = ACLMS(v_in', wind(w, :), stepsize(st), M(i));

        se_CLMS(w,i) = mean(abs(e_CLMS{w,i}).^2);
        se_ACLMS(w,i) = mean(abs(e_ACLMS{w,i}).^2);
    end
end
figure(); 
names = ["High", "Medium", "Low"];
for w = 1:nWinds
    subplot(1, nWinds, w); hold on
    plot(10*log10(se_CLMS(w,:)), 'LineWidth', 1.3);
    plot(10*log10(se_ACLMS(w,:)), 'LineWidth', 1.3);
    legend('CLMS', 'ACLMS');
    title(names(w) + " Wind (\mu = " + stepsize(w) + ")", 'FontSize', 13)
    xlabel('Filter Order', 'FontSize', 13);
    ylabel('Mean Square Error (dBs)', 'FontSize', 13)
end
% overfitting happens more with ACLMS (see larger orders) bc high number of
% degrees of freedom

