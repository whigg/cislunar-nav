function df = grad_central(f, x, h)
%GRAD_CENTRAL computes the gradient of a function, f, at the specified 
%state x using the central difference method.
%   Input:
%    - f; function handle
%    - x; state vector to compute the gradient around
%    - h; step size

n = length(x);
df = zeros(size(x));
for i=1:n
    step = zeros(size(x));
    step(i) = h;
    df(i) = (f(x + step) - f(x - step)) / (2*h);
end
end

