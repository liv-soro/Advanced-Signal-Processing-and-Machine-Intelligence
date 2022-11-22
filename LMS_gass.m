function [weights_gass, prediction, error] = LMS_gass(x, mdl_inputs, order, stepsize, rho, leak, algorithm, alpha)
%LMS_GASS Computes the filter coefficients (weights), prediction and prediction
%error according to the Least Mean Square Algorithm with gradient adaptive
%step size
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
    weights_gass = zeros(nSamp +1, order);
    prediction = zeros(nSamp, 1);
    error = ones(nSamp, 1);
    psi = zeros(nSamp +1, order);
    stepsizes = zeros(nSamp+1, 1);
    stepsizes(1) = stepsize;

    for n = 1:nSamp
        
        % Predicted filter output
        prediction(n, 1) = weights_gass(n, :) * mdl_inputs(n, :)'; % y(n) is the inner product of the delayed input and the current weights
        
        % Prediction Error
        error(n) = x(n) - prediction(n);
        
        % Update predicted weights: w = w + mu*e*x;
        weights_gass(n+1, :) = (1 - stepsizes(n) * leak) * weights_gass(n, :) + stepsizes(n) * error(n) * mdl_inputs(n, :); 
        
        % Update adaptive step size
        stepsizes(n+1) = stepsizes(n) + rho * error(n) * mdl_inputs(n, :) * psi(n, :)'; %scambiato i transpose
        
        switch algorithm
            case 'ben'
                psi(n+1, :) = ((eye(order) - stepsizes(n) * mdl_inputs(n, :)' * mdl_inputs(n, :)) * psi(n, :)' + error(n) * mdl_inputs(n, :)')'; % un casino con i transpose (fai l'opposto (psi = zeros(order, nSamp)))
            case 'af'
                psi(n+1, :) = alpha * psi(n, :) + error(n) * mdl_inputs(n, :);
            case 'mx'
                psi(n+1, :) = error(n) * mdl_inputs(n, :);
            otherwise
                nothing = 0;
        end
    end
end
