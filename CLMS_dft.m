function [weights_h, prediction, error] = CLMS_dft(input, signal, stepsize, leak)
    % DFT as a least square solution using the Complex Least Mean Square 
    % Algorithm 
    % Inputs:
    %   signal: output signal (desired)
    %   input: input signal (row vector)
    %   stepsize: is the initial learning rate for an adaptive algorithm
    %   and the learning rate/step size for a standard algorithm
    %   M: filter length
    %   leak: leakage coefficient
    % Outputs:
    %   weights_h : standard coefficients (on original signal)
    %   prediction: filter output (estimated true data)
    %   error: prediction error
    
    [~, nCols] = size(input);
    if nCols == 1
        input = input';
    end
    

    nSamp = length(input);
    weights_h = zeros(nSamp, nSamp + 1, 'like', 1i);
    error = zeros(1, nSamp, 'like', 1i);
    prediction = zeros(1, nSamp, 'like', 1i);
    
    for n = 1:nSamp
        
        % Predicted Filter Output
        prediction(n) = weights_h(:, n)' * input(:, n);
        
        % Prediction Error
        error(n) = signal(n) - prediction(n);

        % Update Predicted Weights: w = w + mu*e*x;
        weights_h(:, n+1) = (1 - stepsize * leak) * weights_h(:, n) + stepsize * conj(error(n)) * input(:, n); 
        
    end
end