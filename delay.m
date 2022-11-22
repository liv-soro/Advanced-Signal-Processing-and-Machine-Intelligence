function delayed = delay(x, order)
    [nRows, ~] = size(x);
    if nRows == 1
        x = x';
    end
    nSamp = length(x);
    delayed = zeros(nSamp, order);
    
    for o = 1:order
       delayed(:, o) = [zeros(o, 1); x(1:nSamp - o, 1)]; 
    end
end