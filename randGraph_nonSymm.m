function [A,Abinary] = randGraph_nonSymm(n,selfLoops,type,p)
% [A,Abinary] = randGraph_nonSymm(n,selfLoops,type,p)
% function for generating random connected graphs using 
% Erdos-Renyi random graph model
%
% Outputs
% @A = (nonsymmetric) weighted, row-stochastic adjacency matrix A\in[0,1]^{n\times n}
% @Abinary = binary adjacency matrix A\in{0,1}^{n\times n}
%
% Inputs
% @n = dimension of matrix to be generated (number of nodes)
% @selfLoops = {true, false} if true, then A will have strictly positive
%   diagonal entries, otherwise diagonal entries will be 0
% @type = {'complete','cycle','irreducible'} specifies type of graph structure.
%   'full' = complete graph (all entries strictly positive)
%   'cycle' = directed ring graph
%   'irreducible' = [DEFAULT] strongly connected graph (or irreducible 
%       matrix) generated using
%       Erdos-Renyi random graph model with connectivity p
% @p = connectivity probability p, if type = 'primitive', p must be
%   specified

switch type
    case 'complete'
        Abinary = ones(n,n) - eye(n,n);
    case 'cycle' % cycle graph case 
        Abinary = diag(ones(n-1,1),1) + diag(1,-(n-1)); 
    otherwise % type = 'primitive', generate erdos-renyi random graph
        % ensure that the graph is connected: throw out reducible cases
        irreducible = false;
        if p > 1
            error('Input ''p'' must be less than or equal to 1.');
        end
        while irreducible == false 
            Abinary = double(rand(n,n) < p);
            Abinary = Abinary - diag(diag(Abinary));
            sumTemp = zeros(n);
            for k = 0:1:n-1
                sumTemp = sumTemp + Abinary^k;
            end
            irreducible = all(sumTemp,'all'); % check connectivity
        end
end %switch

% remove self loops if no self loops is specified (0s on diagonal)
if selfLoops == 1
    Abinary = Abinary + eye(n);
end

% Assign coupling weight to each edge with aij = aji in [1,10]
A = rand(n,n).*Abinary;
A = diag(A*ones(n,1))\A; % Scale everything to be row-stochastic

end