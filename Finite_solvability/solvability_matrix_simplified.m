%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% Implementation of new simpler method for Finite Solvability

function [A,m]=solvability_matrix_simplified(G,C,use_independent)
% Input: 
% C = camera centres (4 x #cameras)
% G = viewing graph
%
% Output:
% A = finite solvability matrix
% m = number of edges in G

if nargin<3
    use_independent=true;
end

%% Properties of the graph

n=numnodes(G); % number of nodes in input graph
m=numedges(G); % number of edges in input graph
deg=degree(G); % degree of nodes

%% build the REDUCED finite solvability matrix

disp('Building the reduced solvability matrix...')
if use_independent
    nEq = sum(11*(deg-1)); % number of equations
else
    nEq = sum(20*(deg-1)); % number of equations
end

A = spalloc(nEq,16*m,nEq*16*2); % initialization (#rows, #columns, #nonzero)

k=1; % index current equation 
for i=1:n % for all the nodes ... 
    % let us consider node i in the input graph
    
    % build local equations (depend on camera centers)
    if use_independent
        Ai=finite_equations_independent(C(:,i)); % 20 equations (redundant)
    else
        Ai=finite_equations(C(:,i)); % 11 equations (independent)
    end
    
    N = neighbors(G,i);  % find nieghbors of node i 
    assert(length(N)==deg(i)) % sanity check
    
    tau=findedge(G,i,N(1)); % first neighbor is the root
    for j=2:deg(i)
        nu=findedge(G,i,N(j)); % pick another edge in the tree
        
        % update the global matrix
        if use_independent
            A(11*k-10:11*k,16*tau-15:16*tau)=Ai;
            A(11*k-10:11*k,16*nu-15:16*nu)=-Ai;            
        else
            A(20*k-19:20*k,16*tau-15:16*tau)=Ai;
            A(20*k-19:20*k,16*nu-15:16*nu)=-Ai;
        end
        
        k=k+1;
    end
    
end

end
