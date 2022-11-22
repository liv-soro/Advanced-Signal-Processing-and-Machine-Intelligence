function [circ_quotient, circ_coeff] = circularity(signal)
    % CIRCULARITY Computes the circularity quotient and the ciruclarity
    % coefficient of a signal (which quantify the signal's degree of
    % non-circularity
    % Inputs:
    %   signal: complex vector
    % Outputs
    %   circ_quotient: pseudocovariance/covariance -> complex
    %   circ_coeff: abs(pseudocovariance)/covariance -> real

    c = mean(abs(signal).^2); % covariance
    p = mean(signal.^2); % pseudocovariance

    circ_quotient = p/c;
    circ_coeff = abs(p)/c;

end