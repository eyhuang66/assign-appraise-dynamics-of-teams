function [tall,vall,wall,Pall] = get_numericalSoln(A0,w0,Pfunc,tend,cutoff)
% [tall,yall,Pall,Aend,colA_idx] = getSoln(A0,w0,Pfunc,cutoff)
% This function computes the numerical solution for the reduced order
% model, (v(t),w(t))
%
% outputs = [tall,yall,Pall,Aend,colA_idx], 
%   where appraisal matrix at index t is v(t) = y(t,1:n),[n,n])' (column vector)
%   where workload at index t is w(t) = y(t,n+1:end)' (column vector)
%   @tall: time vector
%   @yall: numerical solution of reduced order appraisal states and work
%   @Pall: individual performance function values at every time step
%
% inputs = (A0,w0,Pfunc,cutoff)
%   @A0: nxn initial appraisal matrix
%   @w0: nx1 initial workload assignment
%   @Pfunc: nx1 individual performance functions 
%   @tend: solve numerically for time units [0,tend]  
%   @cutoff: maximum time [seconds] function is allowed to run for, in      
%   case solver takes too long

n = length(A0);

dynamics = @(t,y) model_dc_reducedOrder(t,y,Pfunc,A0);

% initialize initial condition
y0 = [ones(n,1); w0]; 

% setup storage variables for each time step
yall = y0';
Pall = Pfunc(w0)';
tall = 0;

tstep = 0.5;
tset = 0;
err = 10; 
tic
while err > 1e-5
    tset = tset + 1;
    % time interval [0,tend] is split up to improve solver speed
    tspan = [(tset-1)*tstep,tset*tstep]; 

    options = odeset('MaxStep',tstep/2,...
        'RelTol',1e-12,...
        'NonNegative',1:1:2*n);
    [t,y] = ode113(dynamics,tspan,y0,options);

    % populate variables for all time steps
    tall = [tall; t(2:end)];
    yall = [yall; y(2:end,:)];
    Pall = [Pall; Pfunc(y(2:end,end-n+1:end)')'];

    % update new initial values for next tset
    y0 = y(end,:);

    err = norm(yall(end,:)'-yall(end-3,:)',inf);

    time = toc;
    if time > cutoff 
        fprintf('get_numericalSoln has been running for too long\n');
        toc
        err = 1e-20;
    end
    
    if tall(end) > tend 
        err = 1e-20;
    end
end %tset
vall = yall(:,1:n); 
wall = yall(:,n+1:end);
end %function