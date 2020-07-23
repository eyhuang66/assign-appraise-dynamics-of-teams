% example script for running reduced order dynamics
plotFigs = 1; % 1 to plot numerical soln. 0 to supress plot.

% define team parameters
n = 6; % number of nodes
tol = 5e-3; % ensure initial values are strictled bounded from 0

% specify end time for ode45
tend = 1e3;

% monte carlo parameters
randTot = 27000;
% randTot = 10;

% generate random performance functions
skillTest = 0;
while skillTest == 0
    skill = rand(n,1); 
%         skill(skill < tol) = tol;
    skill = skill./sum(skill);
    if min(skill) >= tol
        skillTest = 1;
    end
end
gamma = max(tol,rand(1))*ones(n,1);
Pfunc = @(x) (skill./x).^gamma;

% randomly generate initial appraial matrix
% Erdos-renyi random graph model used to generate the topology
A0 = randGraph_nonSymm(n,true,'irreducible',0.3);
[r1,c1] = find(A0 < tol);
[r2,c2] = find(A0 > 0);
for ii = 1:1:length(r1)
    if ismember([r1(ii),c1(ii)],[r2,c2],'rows') == 1
        A0(r1(ii),c1(ii)) = tol;
    end
end
A0 = diag(A0*ones(n,1))\A0;

% generate random initial work assignment
w0 = rand(n,1); 
w0(w0 < tol) = tol;
w0 = w0./sum(w0);

% numerically solve for solution of system
% cap solver time at 30min = 1800sec
warning('off','all')
[t,v,w,pall] = get_numericalSoln(A0,w0,Pfunc,tend,1800);

if plotFigs == 1
    figure();
    subplot(3,1,1); plot(t,pall); title('Performance, p(w(t))');
    subplot(3,1,2); plot(t,v); title('Reduced order appraisal states, v(t)');
    subplot(3,1,3); plot(t,w); title('Workload, w(t)');
end

% convert reduced order state at final time tend back to full appraisal
% state
fprintf('Appraisal state at tend: \n');
Aend = diag(A0*v(end,:)')\A0*diag(v(end,:))