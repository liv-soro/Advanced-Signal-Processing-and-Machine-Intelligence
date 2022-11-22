function [weights_h, weights_g, prediction, error] = ACLMS(input, signal, stepsize, M)
    % Augmented Complex Least Mean Square Algorithm (accounts for improper
    % data)
    % Inputs:
    %   true_data: output signal
    %   input: input signal (row vector)
    %   stepsize: is the initial learning rate for an adaptive algorithm
    %   and the learning rate/step size for a standard algorithm
    % Outputs:
    %   weights_h : standard coefficients (on original signal)
    %   weights_g : conjugate coefficients 
    %   prediction: filter output (estimated true data)
    %   error: prediction error

    [~, nCols] = size(input);
    if nCols == 1
        input = input';
    end
    nSamp = length(input);
    weights_h = zeros(M, nSamp + 1, 'like', 1i);
    weights_g = zeros(M, nSamp + 1, 'like', 1i);
    error = zeros(1, nSamp, 'like', 1i);
    prediction = zeros(1, nSamp, 'like', 1i);
    
    for n = 1:nSamp
        
        % Predicted Filter Output
        prediction(n) = weights_h(:, n)' * input(:, n) + weights_g(:, n)' * conj(input(:, n));
        
        % Prediction Error
        error(n) = signal(n) - prediction(n);

        % Update Predicted Weights: w = w + mu*e*x;
        weights_h(:, n+1) = weights_h(:, n) + stepsize * conj(error(n)) * input(:, n); 
        weights_g(:, n+1) = weights_g(:, n) + stepsize * conj(error(n)) * conj(input(:, n)); 
    end
end