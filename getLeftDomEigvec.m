function [x] = getLeftDomEigvec(A)
% [x] = getLeftDomEigvec(A)
% get normalized left dominant eigenvector of square matrix A
% sum(x) = 1

if size(A,1) ~= size(A,2)
    error('getLeftDomEigvec(A): A is not a square matrix.\n');
end
[V,D] = eig(A');
[~,idx_spec] = max(abs(diag(D)));
vleft = V(:,idx_spec); 
x = vleft./sum(vleft);
if norm(A'*x-x,inf) > 1e-6
    error('getLeftDomEigvec found incorrect eigenvector');
end %if
end