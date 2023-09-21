%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% TEST for finite solvability

function [issolvable,lambda]=finite_solvability(S,nE,method)
% Input:
% S = solvability matrix
% nE = nodes in the line graph = edges in input graph
% method = 'rank' (full) rank - very slow for large sparse matrices
% method = 'eigs' use eigs
% method = 'svds' use svds
%
% Output:
% issolvable = true if the graph finite solvable, false otherwise
% m = number of edges in G
% lambda = value that when under a threshols determines a rank drop;
%       the meaning depends on the method and can be [];


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

%% Check the rank

% default tol. Formula taken from SuiteSparseQR (similar to Matlab's rank)
% tol = max(m,n)*eps(normest(A,1e-4));

tol = 1e-10;

switch method

    case 'rank'
        disp('Computing the full rank...')
        % size(A)
        r =  rank(full(A));
        issolvable =( size(A,2) == r );
        lambda=[];

    case 'svds'
        disp('Computing the smallest singular value with svds...')
        lambda = svds(A,1,'smallest','Display',true); % based on QR
        %lambda = svds(A,1,tol,'Display',true); % uses EIGS of matrix [0 A; A' 0] --> SLOW
        issolvable=lambda(end)>tol;

    case 'eigs'
        disp('Computing the smallest eigenvalue with eigs...')
        M=A'*A;
        assert(issparse(M));
        %disp(['amount of entries = ' num2str( nnz(M)/(size(M,1)*size(M,2))*100 ) '%']);
        lambda = eigs(M,1,tol,'Display',true); % works for non-solvable
        % lambda = eigs(M,1,'sm','Display',true); this explodes for non-solvable
        issolvable=lambda(end)>tol;

    otherwise
        error('method not implemented')
end


end


