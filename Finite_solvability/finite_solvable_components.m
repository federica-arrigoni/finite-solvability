%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% COMPUTATION of maximal finite-solvable components

function [fcomp,U]=finite_solvable_components(S,nE,method)
% Input:
% S = solvability matrix
% nE = nodes in the line graph = edges in input graph
% method = 'svd' use Matlab function null (based on svd) --> small graphs
% method = 'nulls' use  Gotsman-Toledo nulls --> large sparse graphs
%
% Output:
% fcomp = index of solvable components (vector of length #edges)
% U = null-space of S


%% Fix projective & scale ambiguity

disp('Fixing the trivial d.o.f.')

% fix first transformation to be something
B=sparse(16,16*nE);
B(1:16,1:16)=speye(16);

% fix the scale of remaining transformations (entries sum to 1)
%     D=sparse(nL-1,16*nL);
%     for k=2:nL
%         D(k-1,16*k-15:16*k)=1;
%     end

% faster version
i     = repmat(1:(nE-1),16,1);
D     = sparse(i(:), 17:nE*16, ones(1,16*(nE-1)));

A=[S;B;D];
clear S B D;

%% Compute the null-space

% default tol. Formula taken from SuiteSparseQR (similar to Matlab's rank)
% tol = max(m,n)*eps(normest(A,1e-4));

tol = 1e-10;

[m,n] = size(A);
assert(m>=n); % Not all the methods work with slim matrices

%% Compute the null-space
switch method

    case 'svd' % full matrix
        disp('Computing null-space with null(full(A))...')
        U = null(full(A));
        dim_null =  size(U,2);

    case 'nulls' % Best so far
        disp('Computing null-space with Gotsman-Toledo...')
        [U,flag,dim_null] = gt_nulls(A,1,tol);

        if flag==1
            warning('the null space might be larger' )
        end

    otherwise
        error('method not implemented')

end

assert(dim_null ==  size(U,2)) ;  % dimension of null space

%% Find all components: cluster rows of U to assign components

tol=1e-4; % to decide if two columns are the same

if size(U,2) == 0
    % Solvable, only ONE component
    fcomp =  ones(nE,1);
else
    % This is U where every 16 rows are summed up to one single row (abs value)
    collapsedU = reshape(accumarray(kron([1:nE*size(U,2)]', ones(16,1)),abs(U(:))), nE, []);

    [fcomp,C] = myclustering(collapsedU, tol);
    assert(size(C,1) == max(fcomp)); % # of components

end

end


