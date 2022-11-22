function [weights, prediction, error] = LMS_leaky(x, order, stepsize, leakage)
%LMS Computes the 
    % Inputs:
    %   x: Input vector of the process to estimate (M x 1)
    %   order: Order of the model - determines the number of parameters to be
    %   estimated
    %   stepsize: step size (= learning rate) of the adaptive algorithm
    %   leakage: leakage coefficient (0 < leakage < 1) 
    % Outputs: 
    %   weights: Series of weights of the LMS algorithm
    %   prediction: Filter output
    %   error: Prediction error
    
    nSamp = length(x);
    weights = zeros(nSamp +1, order);
    prediction = zeros(nSamp, 1);
    error = zeros(nSamp, 1);
    x_delay = delay(x, order);

    for n = 1:nSamp
        
        % Predicted filter output
        prediction(n, 1) = weights(n, :) * x_delay(n, :)'; % y(n) is the inner product of the delayed input and the current weights
        
        % Prediction Error
        error(n) = x(n) - prediction(n);
        
        % Update predicted weights: w = w + mu*e*x;
        weights(n+1, :) = (1 - stepsize * leakage) * weights(n, :) + stepsize * error(n) * x_delay(n, :); 
        
    end
end

