function estimate = ACF(x, k, type)
% x is the data
N = length(x);
estimate = 0;
for n = k+1:N
    estimate = estimate + x(n)*conj(x(n-k));
end

if strcmpi(type, "unbiased")
    estimate = estimate/(N-k);
elseif strcmpi(type, "biased")
    estimate = estimate/N;
else
    print('prob in ACF function')
end
end