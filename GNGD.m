function [weights_GNGD, prediction, error] = GNGD(x, mdl_inputs, order, stepsize, rho)
%GNGD Computes the filter coefficients (weights), prediction and prediction
%error according to the Generalised Normalised Gradient Descent algorithm
%(adaptive Normalised Least Mean Square Algorithm)
    % Inputs:
    %   x: Input vector of the process to estimate (M x 1)
    %   order: Order of the model - determines the number of parameters to be
    %   estimated
    %   stepsize: step size (= learning rate) of the adaptive algorithm
    % Outputs: 
    %   weights: Series of weights of the LMS algorithm
    %   prediction: Filter output
    %   error: Prediction error
    
    nSamp = length(x);
    weights_GNGD = zeros(nSamp +1, order);
    prediction = zeros(nSamp, 1);
    error = ones(nSamp, 1);
    epsi = zeros(nSamp + 1, 1);
    epsi(1) = 1/stepsize;

    for n = 1:nSamp
        
        % Predicted filter output
        prediction(n, 1) = weights_GNGD(n, :) * mdl_inputs(n, :)'; % y(n) is the inner product of the delayed input and the current weights
        
        % Prediction Error
        error(n) = x(n) - prediction(n);
        
        % Update predicted weights: w = w + mu*e*x;
        weights_GNGD(n+1, :) = weights_GNGD(n, :) + (1/(epsi(n) + norm(mdl_inputs(n,:)).^2)) * error(n) * mdl_inputs(n, :); 
        
        % Update adaptive regularisation factor
        if n>1
            num = error(n)*error(n-1)*dot(mdl_inputs(n,:), mdl_inputs(n-1, :));
            den = (epsi(n-1) + norm(mdl_inputs(n-1,:)).^2).^2;
            epsi(n+1) = epsi(n) - (rho*stepsize) * (num/den);
        end
    end
end