
%% COMPUTATION of maximal finite-solvable components
%% Iteratively identify all nodes that are in the same component as the two 
%% ones used to fix the global projective ambiguity

%% NOTE: this method is highly INEFFICIENT at it re-computes the 
%% solvability matrix for each component. Instead, we should keep the same 
%% solvability matrix and just updates the blocks related to ambiguities 
%% (we should keep track of cameras/fundamental matrices)

function [fcomp,U]=direct_solvable_components(S,G,nE,nV,method,node1,node2)
% Input:
% S = solvability matrix
% nE = edges in input graph
% method = 'svd' use Matlab function null (based on svd) --> small graphs
% method = 'nulls' use  Gotsman-Toledo nulls --> large sparse graphs
%
% Output:
% fcomp = index of solvable components (vector of length #edges)
% U = null-space of S

%% Parameters & Initialization 

% tol = max(m,n)*eps(normest(A,1e-4)); % default tol. Formula taken from SuiteSparseQR (similar to Matlab's rank)
tol_null = 1e-10; % nulls

tol_rows=1e-4; % to decide if two rows are the same

fcomp = nan(nE,1); % initialize components
E=G.Edges{:,1}; % all the edges
ind_comp = 1; % index of current component

%% Find the first component

% Compute the null-space
[U]=nullspace_components(S,method,tol_null);

% Identify all nodes that are in the same component as the two ones used
% to fix the global projective ambiguity
[fnodes]=find_component_2nodes(U,nV,tol_rows,node1,node2);

% Find corresponding edges
fedges=ismember(E,fnodes);
fedges=fedges(:,1) & fedges(:,2); % edges having both endpoints in the component
fcomp(fedges)=ind_comp; % update CURRENT COMPONENT

%% Find other components

while(~isempty(find(isnan(fcomp))))

    % Increment the number of components
    ind_comp=ind_comp+1;

    % find one edge not assigned to the previous components
    ind_edge=find(isnan(fcomp)); ind_edge=ind_edge(1);
    node1=E(ind_edge,1); node2=E(ind_edge,2);

    % recompute the solvability matrix using such nodes to fix ambiguity
    [S]=solvability_matrix_skew(G,true,true,true,node1,node2);

    % Compute the null-space
    [U]=nullspace_components(S,method,tol_null);

    % Identify all nodes that are in the same component as the two ones used
    % to fix the global projective ambiguity
    [fnodes]=find_component_2nodes(U,nV,tol_rows,node1,node2);

    % Find corresponding edges
    fedges=ismember(E,fnodes);
    fedges=fedges(:,1) & fedges(:,2);
    fcomp(fedges)=ind_comp; % CURRENT COMPONENT

end

end



%% Function to compute nullspace of solvability matrix

function [U,dim_null]=nullspace_components(S,method,tol)

[m,n] = size(S);
assert(m>=n); % Not all the methods work with slim matrices

switch method

    case 'svd' % full matrix
        disp('Computing null-space with null(full(A))...')
        U = null(full(S));
        dim_null =  size(U,2);

    case 'nulls' % Best so far
        disp('Computing null-space with Gotsman-Toledo...')
        [U,flag,dim_null] = gt_nulls(S,1,tol);

        if flag==1
            warning('the null space might be larger' )
        end

    otherwise
        error('method not implemented')

end

assert(dim_null ==  size(U,2)) ;  % dimension of null space

end

%% function to identify zero rows (up to threshold)

function [fnodes]=find_component_2nodes(U,nV,tol_rows,node1,node2)

if size(U,2) == 0
    % Solvable, only ONE component
    % fcomp =  ones(nE,1);
    fnodes=1:nV;
else

    U=U(1:12*nV,:); % remove parts related to auxiliary variables
    % TO DO: use auxiliary variables to assert that everything is fine

    % This is U where every 12 rows are summed up to one single row (abs value)
    collapsedU = reshape(accumarray(kron([1:nV*size(U,2)]', ones(12,1)),abs(U(:))), nV, []);

    % Find nodes with zero rows
    fnodes=find(vecnorm(collapsedU, 2, 2) < tol_rows );    
    % Check that the two nodes used to fix ambiguity actually belong to this component
    assert(sum(ismember(fnodes,node1))==1)
    assert(sum(ismember(fnodes,node2))==1)
    
end

end

