
%% COMPUTATION of maximal finite-solvable components
%% Iteratively identify all nodes that are in the same component as the two
%% ones used to fix the global projective ambiguity

%% NOTE: this method is highly INEFFICIENT at it re-computes the
%% solvability matrix for each component. Instead, we should keep the same
%% solvability matrix and just updates the block related to ambiguities

function [fcomp,U]=fast_direct_solvable_components(S,G,nE,nV,method,node1,node2,verbose)
% Input:
% S = solvability matrix
% nE = edges in input graph
% method = 'svd' use Matlab function null (based on svd) --> small graphs
% method = 'nulls' use  Gotsman-Toledo nulls --> large sparse graphs
%
% Output:
% fcomp = index of solvable components (vector of length #edges)
% U = null-space of S

if nargin<8
    verbose=false;
end

%% Parameters & Initialization

% tol = max(m,n)*eps(normest(A,1e-4)); % default tol. Formula taken from SuiteSparseQR (similar to Matlab's rank)
tol_null = 1e-10; % nulls
tol_rows=1e-4; % to decide if two rows are the same

fcomp = nan(nE,1); % initialize components
E=G.Edges{:,1}; % all the edges
ind_comp = 0; % index of current component

%% Find the first component

% Compute the null-space
[U]=nullspace_components(S,method,tol_null);

% Identify all nodes that are in the same component as the two ones used
% to fix the global projective ambiguity
[fnodes]=find_component_2nodes(U,nV,tol_rows,node1,node2);

% Find corresponding edges
[fcomp,ind_comp]=update_component(E,fnodes,fcomp,ind_comp);

%% Find other components

while(~isempty(find(isnan(fcomp))))

    % remove edges already assigned
    H = rmedge(G,find(~isnan(fcomp)));
    deg=degree(H); nodes_remaining=find(deg>0);
    H = subgraph(H,nodes_remaining);

    % find connected components in H
    bins = conncomp(H); ncomp=max(bins);

    % consider each connected component of the subgraph separately
    for k=1:ncomp

        if verbose
            disp(['comp ' num2str(k)])
        end

        nodesH_comp=find(bins==k); % nodes in H in current connected component
        Hconn = subgraph(H,nodesH_comp); % subgraph of H (which in turn is subgraph of G)

        % find solvable components within the current connected component
        [fcomp,ind_comp]=finite_connected_component(Hconn,E,nodes_remaining(nodesH_comp),fcomp,ind_comp,method,tol_null,tol_rows,verbose);

    end
end

end
%%

%% Find solvable components within a connected subgraph of the original graph
%% H is the subgraph, E are edges in original graph 
%% nodes_remaining are indices of nodes in the subgraph wrt the original graph
function [fcomp,ind_comp]=finite_connected_component(H,E,nodes_remaining,fcomp,ind_comp,method,tol_null,tol_rows,verbose)

%figure, plot(H)
nV_H=numnodes(H);
nE_H=numedges(H);

if nE_H==1 % check if the subgraph is made of a single edge

    fnodes=nodes_remaining; % Go back to the original graph
    assert(length(fnodes)==2) % single edge
    [fcomp,ind_comp]=update_component(E,fnodes,fcomp,ind_comp);

    if verbose
        figure, plot(H)
        disp('edge'); fnodes
    end

elseif nE_H==nV_H-1 && ~hascycles(H) % check if the subgraph is a spanning tree

    % edges in the tree: each edge is a component
    fnodesH=H.Edges{:,1};
    fnodes=nodes_remaining(fnodesH); % Go back to the original graph
    [fcomp,ind_comp]=update_component(E,fnodes,fcomp,ind_comp);

    if verbose
        figure, plot(H)
        disp('edge'); fnodes
    end


else
    % select the nodes with highest degree to fix projective ambiguity
    deg=degree(H); [~,node1]=max(deg);
    N = neighbors(H,node1); [~,ii]=max(deg(N)); node2=N(ii);

    % Compute the solvability matrix using such nodes to fix ambiguity
    % working on the subgraph
    [S]=solvability_matrix_skew(H,true,true,true,node1,node2);

    % Compute the null-space
    [U]=nullspace_components(S,method,tol_null);

    % Identify all nodes that are in the same component as the two ones used
    % to fix the global projective ambiguity
    [fnodesH]=find_component_2nodes(U,nV_H,tol_rows,node1,node2);

    % Go back to the original graph
    fnodes=nodes_remaining(fnodesH);
    [fcomp,ind_comp]=update_component(E,fnodes,fcomp,ind_comp);

    if verbose
        figure, plot(H)
        disp('edge'); fnodes
    end


end


end

%% Function to update indices of solvable components

function [fcomp,ind_comp]=update_component(E,fnodes,fcomp,ind_comp)

if min(size(fnodes))==1 % single component

    % Find corresponding edges
    fedges=ismember(E,fnodes);
    fedges=fedges(:,1) & fedges(:,2);

    % Increment the number of components
    ind_comp=ind_comp+1;
    fcomp(fedges)=ind_comp; % CURRENT COMPONENT

else % many components (e.g., many edges in a tree, each one is a component)

    for k=1:size(fnodes,1)
        current_nodes=fnodes(k,:);

        % Find corresponding edges
        fedges=ismember(E,current_nodes);
        fedges=fedges(:,1) & fedges(:,2);

        % Increment the number of components
        ind_comp=ind_comp+1;
        fcomp(fedges)=ind_comp; % CURRENT COMPONENT

    end
end

end


%% Function to compute nullspace of solvability matrix

function [U,dim_null]=nullspace_components(S,method,tol)

switch method

    case 'svd' % full matrix
        disp('Computing null-space with null(full(A))...')
        U = null(full(S));
        dim_null =  size(U,2);

    case 'nulls' % Best so far

        [m,n] = size(S);
        assert(m>=n); % Does not work with slim matrices
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

