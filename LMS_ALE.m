function [weights, prediction, error] = LMS_ALE(signal, stepsize, delta, M)
%LMS Computes the filter coefficients (weights), prediction and prediction
%error using the Adaptive Line Enhancer Algorithm
    % Inputs:
    %   x: Input vector of the process to estimate 
    %   stepsize: step size (= learning rate) of the adaptive algorithm
    %   delta: smallest delay
    %   M: ALE filter length
    % Outputs: 
    %   weights: Series of weights of the LMS algorithm
    %   prediction: Filter output
    %   error: Prediction error
    
    nSamp = length(signal);
    weights = zeros(nSamp +1, M);
    prediction = zeros(size(signal));
    error = zeros(size(signal));
    
    for n = 1:nSamp
        
        % Delay signal
        u = delay2(signal, delta, M);
        
        % Predicted filter output
        prediction(n, 1) = weights(n, :) * u(n, :)'; % y(n) is the inner product of the delayed input and the current weights
        
        % Prediction Error
        error(n) = signal(n) - prediction(n);
        
        % Update predicted weights: w = w + mu*e*x;
        weights(n+1, :) = weights(n, :) + stepsize * error(n) * u(n, :); 
        
    end
end
function delayed = delay2(in, delta, M)
    % del is delta
    % M: filter order
    nSamp = length(in);
    delayed = zeros(nSamp, M);
    
    for k = 1:M
       delay = delta + k - 1;
       delayed(:, k) = [zeros(delay, 1); in(1:nSamp - delay, 1)]; 
    end
end
