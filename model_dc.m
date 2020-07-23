function dydt = model_dc(t,y,Pfunc,A0)
% dydt = model_reducedOrder(t,y,Pfunc)
% function handle for ODE of the full order ASAP model with
% donor-controlled work flow
% dy/dt = f(t,y)
%
% Output
% @dydt = dt/dt, 
% where y = [x;w] = [columns of appraisal matrix stacked; workload]
%
% Input
% @t = time variable (unused since time-invariant system)
% @y = state variable with y = [v; w]
% @Pfunc = performance functions used in the model
% @A0 = initial appraisal matrix

n = size(A0,1);

% A dynamics
A = reshape(y(1:n*n),[n,n]); 
w = y(n*n+1:end);
F = Pfunc(w);
dAdt = A.*( ones(n,1)*F' - (A*F)*ones(1,n) );

% x dynamics
dwdt = -(eye(n)-A')*w;

% reshape back into column vector
dydt = [reshape(dAdt,[n*n,1]); dwdt];
end