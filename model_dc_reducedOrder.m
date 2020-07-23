function dydt = model_dc_reducedOrder(t,y,Pfunc,A0)
% dydt = model_dc_reducedOrder(t,y,Pfunc,A0)
% function handle for ODE of the reduced order ASAP model with
% donor-controlled work flow
% dy/dt = f(t,y)
%
% Output
% @dydt = dy/dt, where y = [v;w] = [reduced appraisal; workload]
%
% Input
% @t = time variable (unused since time-invariant system)
% @y = state variable with y = [v; w]
% @Pfunc = performance functions used in the model
% @A0 = initial appraisal matrix

n = size(A0,1);
v = y(1:n); w = y(end-n+1:end);
A = diag(A0*v)\A0*diag(v); % map v coordinates back to A states

p = Pfunc(w);

dvdt = v.*(p - (p'*A'*w)*ones(n,1));
dxdt = -w+A'*w;
dydt = [dvdt; dxdt];
end