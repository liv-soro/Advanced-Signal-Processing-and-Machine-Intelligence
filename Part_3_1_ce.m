% 3. Widely Liner Filtering and Adaptive Spectrum Estimation
clear; close all; clc

%% 3.1. Complex LMS and Widely Linear Modelling
% c) Frequency Estimation in Three Phase Power
nPhases = 3;
fs = 3000; % sampling frequency
fo = 50; % nominal frequency
nSamp = 1000; 
phaseshift = [0; - 2*pi/3; 2*pi/3];
phasedist = [0; 0; 0];
clarke = [sqrt(2)/2,    sqrt(2)/2,  sqrt(2)/2; ...
          1,            -1/2,       -1/2; ...
          0,            sqrt(3)/2,  -sqrt(3)/2];
clarke = sqrt(2/3).* clarke;

% Balanced System
n = 0:nSamp-1;
amps = ones(nPhases, 1);
bal_abc = amps .*  cos(2*pi*(fo/fs)*n + phaseshift  + phasedist);
bal_alphabet = clarke * bal_abc;
bal_clarke = bal_alphabet(2, :) + 1i * bal_alphabet(3,:);
bal_clarke_exp = sqrt(3/2) * amps(1) * exp(1i * (2*pi*(fo/fs)*n));

[~, bal_circ] = circularity(bal_clarke);
figure(); subplot(1, 3, 1)
scatter(real(bal_clarke(:)), imag(bal_clarke(:)))
xlabel('Real Part'); ylabel('Imaginary Part');
legend( "|\rho| = " + round(bal_circ*100)/100);
title('Balanced System')
set(gca, 'FontSize', 13)
axis equal
[minlim, maxlim] = set_axis_equal(bal_clarke);
xlim([minlim, maxlim])
ylim([minlim, maxlim])

% Unbalanced System - different constant amplitudes
AMPS_ub = [0.7, 1, 1.2; 0.4, 1, 1.5; 0.1, 1, 1.8]';

subplot(1, 3, 2); hold on
for i=1:3
    amps_ub = AMPS_ub(:,i);
    unbal_abc = amps_ub .*  cos(2*pi*(fo/fs)*n + phaseshift  + phasedist);
    unbal_alphabet = clarke * unbal_abc;
    A = (sqrt(6)/6) * dot(amps_ub, exp(1i * phasedist));
    B = (sqrt(6)/6) * dot(amps_ub, exp(- 1i * (phasedist - phaseshift)));
    unbal_clarke = A * exp(1i *(2*pi*(fo/fs)*n)) + B * exp(-1i *(2*pi*(fo/fs)*n));

    [~, unbal_circ] = circularity(unbal_clarke);
    leg{i} = "|\rho| = " + round(unbal_circ*100)/100;
    scatter(real(unbal_clarke(:)), imag(unbal_clarke(:)))
end
xlabel('Real Part'); ylabel('Imaginary Part');
title({'Unbalanced System', '(different amplitudes)'})
legend(leg); set(gca, 'FontSize', 13)
axis equal
[minlim, maxlim] = set_axis_equal(unbal_clarke);
xlim([minlim, maxlim])
ylim([minlim, maxlim])


% Unbalanced System - different phase distrotions
PHASEDIST_ub = [0, -pi/2, -pi/2; 0, -sqrt(2)*pi/2, sqrt(2)/2; 0, pi, -pi]';

subplot(1, 3, 3); hold on
for i=1:3
    phasedist_ub = PHASEDIST_ub(:,i);
    unbal2_abc = amps .*  cos(2*pi*(fo/fs)*n + phaseshift  + phasedist);
    unbal2_alphabet = clarke * unbal2_abc;
    A = (sqrt(6)/6) * dot(amps, exp(1i * phasedist_ub));
    B = (sqrt(6)/6) * dot(amps, exp(- 1i * (phasedist_ub - phaseshift)));
    unbal2_clarke = A * exp(1i *(2*pi*(fo/fs)*n)) + B * exp(-1i *(2*pi*(fo/fs)*n));

    [~, unbal2_circ] = circularity(unbal2_clarke);
    leg{i} = "|\rho| = " + round(unbal2_circ*100)/100;
    scatter(real(unbal2_clarke(:)), imag(unbal2_clarke(:)))
end
xlabel('Real Part'); ylabel('Imaginary Part');
title({'Unbalanced System', '(different phase distortions)'})
legend(leg)
set(gca, 'FontSize', 13)
axis equal
[minlim, maxlim] = set_axis_equal(unbal2_clarke);
xlim([minlim, maxlim])
ylim([minlim, maxlim])

%% e) Frequency estimation of three-phase voltages pwr system with CLMS and ACLMS
clear; close all; clc
%% 
fs = 1000; % sampling frequency
fo = 50; % nominal frequency
nSamp = 1000;
n = 0:nSamp-1;
nPhases = 3;
phi = 0;
phaseshift = [0; - 2*pi/3; 2*pi/3];

% Balanced System
amps = [1; 1; 1];
v_bal = sqrt(3/2) * amps(1) * exp(1i * (2*pi*(fo/fs)*n + phi));
[~, circ_bal] = circularity(v_bal);

% Unbalanced System 1
amps_ub = [0.3; 0.6; 0.9];
phasedist = [0;0;0];
A = (sqrt(6)/6) * dot(amps_ub, exp(1i * phasedist));
B = (sqrt(6)/6) * dot(amps_ub, exp(- 1i * (phasedist - phaseshift)));
v_unbal1 = A * exp(1i *(2*pi*(fo/fs)*n + phi)) + B * exp(-1i *(2*pi*(fo/fs)*n + phi));
[~, circ_unbal1] = circularity(v_unbal1);

% Unbalanced System 2
phasedist_ub = [0; 0.2*pi; 0.4*pi];
A = (sqrt(6)/6) * dot(amps, exp(1i * phasedist_ub));
B = (sqrt(6)/6) * dot(amps, exp(- 1i * (phasedist_ub - phaseshift)));
v_unbal2 = A * exp(1i *(2*pi*(fo/fs)*n + phi)) + B * exp(-1i *(2*pi*(fo/fs)*n + phi));
[~, circ_unbal2] = circularity(v_unbal2);


%% Estimation of systems using CLMS and ACLMS
st = 0.05; %[1, 0.1, 0.01, 0.001, 0.0001]
bal_in = delay(v_bal', 1);
[h_clms_bal, ~, e_clms_bal] = CLMS(bal_in', v_bal, st, 1);
[h_aclms_bal, g_aclms_bal, ~, e_bal] = ACLMS(bal_in', v_bal, st, 1);

unbal1_in = delay(v_unbal1', 1);
[h_clms_unbal1, ~, e_clms_unbal1] = CLMS(unbal1_in', v_unbal1, st, 1);
[h_aclms_unbal1, g_aclms_unbal1, ~, e_unbal1] = ACLMS(unbal1_in', v_unbal1, st, 1);

unbal2_in = delay(v_unbal2', 1);
[h_clms_unbal2, ~, e_clms_unbal2] = CLMS(unbal2_in', v_unbal2, st, 1);
[h_aclms_unbal2, g_aclms_unbal2, ~, e_unbal2] = ACLMS(unbal2_in', v_unbal2, st, 1);

%%
% Estimate System Frequency
f_clms_bal = f_SL(h_clms_bal, fs);
f_aclms_bal = f_WL(h_aclms_bal, g_aclms_bal, fs);
f_clms_unbal1 = f_SL(h_clms_unbal1, fs);
f_aclms_unbal1 = f_WL(h_aclms_unbal1, g_aclms_unbal1, fs);
f_clms_unbal2 = f_SL(h_clms_unbal2, fs);
f_aclms_unbal2 = f_WL(h_aclms_unbal2, g_aclms_unbal2, fs);


%% Plot
figure(); subplot(3, 1, 1); hold on
plot(abs(f_clms_bal), 'LineWidth', 1.3);
plot(abs(f_aclms_bal), 'LineWidth', 1.3)
plot(0:length(f_clms_bal)-1, fo*ones(1, length(f_clms_bal)), '--', 'LineWidth', 1.3)
legend('CLMS', 'ACLMS', 'True f_o')
title('Frequency Estimation for Balanced System', 'FontSize', 13)
xlabel('Time Step', 'FontSize', 13); ylabel('Frequency (Hz)', 'FontSize', 13); 
ylim([0, 80]); xlim([0, 1000]);

subplot(3, 1, 2); hold on
plot(abs(f_clms_unbal1), 'LineWidth', 1.3);
plot(abs(f_aclms_unbal1), 'LineWidth', 1.3);
plot(0:length(f_clms_bal)-1, fo*ones(1, length(f_clms_bal)), '--', 'LineWidth', 1.3)
legend('CLMS', 'ACLMS', 'True f_o')
title('Frequency Estimation for Unbalanced System (Amplitude)', 'FontSize', 13)
xlabel('Time Step', 'FontSize', 13); ylabel('Frequency (Hz)', 'FontSize', 13); 
ylim([0, 80]); xlim([0, 1000]);

subplot(3, 1, 3); hold on
plot(abs(f_clms_unbal2), 'LineWidth', 1.3);
plot(abs(f_aclms_unbal2), 'LineWidth', 1.3);
plot(0:length(f_clms_bal)-1, fo*ones(1, length(f_clms_bal)), '--', 'LineWidth', 1.3)
legend('CLMS', 'ACLMS', 'True f_o')
title('Frequency Estimation for Unbalanced System (Phase)', 'FontSize', 13)
xlabel('Time Step', 'FontSize', 13); ylabel('Frequency (Hz)', 'FontSize', 13); 
ylim([0, 80]); xlim([0, 1000]);

%% Convergence
disp(f_clms_bal(end))
disp(f_aclms_bal(end))
disp(f_clms_unbal1(end))
disp(f_aclms_unbal1(end))
disp(f_clms_unbal2(end))
disp(f_aclms_unbal2(end))

% st = 1
%   -50.0000
%    NaN
%   -60.1404
%    50.0000
%   -50.0000
%    NaN

% st = 0.1
%   -50.0000
%    50.0000
%   -42.7061
%    50.0000
%   -50.0000
%    50.0000

% st = 0.05
%  -50.0000
%    50.0000
%   -42.6446
%    50.0000
%   -50.0000
%    50.0000

% st = 0.01
%   -50.0000
%    50.0000
%   -42.6810
%    49.9108
%   -50.0000
%    50.0000

% st = 0.001
%   -50.0000
%    50.1586
%   -42.6831
%      0
%   -50.0000
%    50.1586

% st = 0.0001
%   -50.0000
%    50.0339
%   -42.6799
%      0
%   -50.0000
%    50.0339
%%

function fo = f_SL(h, fs)
    fo = (fs/(2*pi)) * atan2(imag(h), real(h));
end

function fo = f_WL(h, g, fs)
    fo = (fs/(2*pi)) * atan2(real(sqrt(imag(h).^2 - abs(g).^2)), real(h));
end






