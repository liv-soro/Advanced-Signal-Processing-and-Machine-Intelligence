function [weights, prediction, error] = LMS_tanh(signal, input, order, stepsize, a, bias, train)
    nSamp = length(signal);
    prediction  = zeros(1, nSamp);
    error  = zeros(1, nSamp);
    weights  = zeros(order+bias, nSamp);
    
    if train
        [weights0,a] = pretrain(signal, input, 100, 20, order, stepsize, a);
        weights  = weights0 .* ones(order+1, nSamp);
        stepsize = 2*10^(-7);
    end

    for i = order:nSamp-1
        if bias
            x = [1; input(i:-1:i-order+1)];
        else 
            x = input(i:-1:i-order+1);
        end
        prediction(i) = a*tanh(weights(:,i).' * x);
        error(i) = signal(i) - prediction(i);
        weights(:,i+1) = weights(:,i) + stepsize * error(i) * x;
    end
end