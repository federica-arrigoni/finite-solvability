%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% Implementation of Trager et al.'s paper (ECCV 2018)

function [A,nL,mL]=solvability_matrix(G,C)
% Input: 
% C = camera centres (4 x #cameras)
% G = viewing graph
%
% Output:
% A = finite solvability matrix
% nL = number of nodes in the line graph = number of edges in G
% mL = number of edges in the line graph (related to number of equations)

%% build the line graph

B=incidence(G); % incidence matrix
AL=spones(B'*B); AL=triu(AL,1); AL=AL+AL'; % adjacency matrix of the line graph
L=graph(AL); % line graph
EL=L.Edges{:,1}; % edges in the line graph
mL=numedges(L); % number of edges in line graph
nL=numnodes(L); % number of nodes in line graph

%% build the finite solvability matrix

disp('Building the solvability matrix...')

A = spalloc(20*mL,16*nL,640*mL); % initialization (#rows, #columns, #nonzero)
for k=1:mL
    tau=EL(k,1); nu=EL(k,2); % edge in the line graph (pair of nodes)
    
    [s1,t1] = findedge(G,tau); % node in the line graph is an edge in the original graph
    [s2,t2] = findedge(G,nu); % node in the line graph is an edge in the original graph
    
    i=intersect([s1,t1],[s2,t2]); % common node (in the original graph)
    
    Ai=finite_equations(C(:,i)); % build local equations
    
    % update the global matrix
    A(20*k-19:20*k,16*tau-15:16*tau)=Ai;
    A(20*k-19:20*k,16*nu-15:16*nu)=-Ai;
end


end





