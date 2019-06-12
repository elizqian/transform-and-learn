% solves ODE with residual f and initial condition x0 using forward Euler

% INPUTS
% f         dy/dt = f
% y0        initial state
% t         vector of times at which to output y

% OUTPUT
% y_all     state data at times specified in t

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

function y_all = forwardEuler(f,y0,t)
K = length(t)-1;
n = length(y0);

y_all = zeros(n,K);
y_all(:,1) = y0 + (t(2)-t(1))*f(y0);
for i = 2:K
    y_all(:,i) = y_all(:,i-1) + (t(i+1)-t(i))*f(y_all(:,i-1));
end