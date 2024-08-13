%% DIRECT SOLVABILITY
%% TEST for finite solvability
%% We assume the projective and scale ambiguity have been fixed and are 
%% already included in the solvability matrix

function [issolvable,lambda]=direct_finite_solvability(A,method)
% Input:
% S = solvability matrix
% nE = edges in input graph
% method = 'rank' (full) rank - very slow for large sparse matrices
% method = 'eigs' use eigs
% method = 'svds' use svds
%
% Output:
% issolvable = true if the graph finite solvable, false otherwise
% m = number of edges in G
% lambda = value that when under a threshols determines a rank drop;
%       the meaning depends on the method and can be [];

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
        disp(['amount of entries = ' num2str( nnz(M)/(size(M,1)*size(M,2))*100 ) '%']);
        lambda = eigs(M,1,tol,'Display',true); % works for non-solvable
        issolvable=lambda(end)>tol;

    case 'eigs_sm'
        disp('Computing the smallest eigenvalue with eigs...')
        M=A'*A;
        assert(issparse(M));
        disp(['amount of entries = ' num2str( nnz(M)/(size(M,1)*size(M,2))*100 ) '%']);
        lambda = eigs(M,1,'sm','Display',true); % this explodes for non-solvable on ICCV23
        issolvable=lambda(end)>tol;

    otherwise
        error('method not implemented')
end


end


