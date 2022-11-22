function [weights, a] = pretrain(signal, input, epochs, samples, order, stepsize, a)
    weights  = zeros(order+1, 1);
    
    for n = 1:epochs
        for i = order:samples-1
            x = [1; input(i:-1:i-order+1)];
            prediction = a*tanh(weights.' * x);
            error = signal(i) - prediction;
            weights = weights + stepsize * error * x;
        end
    end
end