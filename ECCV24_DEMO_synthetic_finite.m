
%% DEMO to test finite solvability on small synthetic cases

%% Adapted from the code of the following paper:
%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023

clc, clear, close all
addpath(genpath('./'))

%% Select one example

% % Graph with 9 nodes that is finite solvable but not solvable (2 solutions)
% load('./GraphsSynthetic/G9.mat'); n=size(A,1);

% % Minimal solvable graph with 90 nodes
% load('./GraphsSynthetic/G90.mat'); n=size(A,1); 

% Random graph
n=20; % nodes
A=rand(n)<0.2;
A=triu(A,1); A=A+A';

%% Plot graph

G=graph(A);
figure, plot(G)

%% Check finite solvability

cc=conncomp(G); % compute connected components

if max(cc)~=1
    disp('Not connected, hence not solvable');
    is_solvable=false;

else

    method='rank'; % FOR LARGE-SCALE EXAMPLES PLEASE USE 'eigs' 
    C=rand(4,n)*100; % random camera centres

    % Trager et al. ECCV 2018
    [S_Trager,nL,mL]=solvability_matrix(G,C);
    [issolvable_Trager,lambda_Trager]=finite_solvability(S_Trager,nL,method);
    disp(['Trager et al.: Is the graph finite-solvable? ' num2str(issolvable_Trager)])
    disp(' ')

    % ICCV 2023
    [S_reduced,nE]=solvability_matrix_simplified(G,C,true);
    [issolvable,lambda]=finite_solvability(S_reduced,nE,method);
    disp(['ICCV 2023: Is the graph finite-solvable? ' num2str(issolvable)])
    disp(' ')

    % Direct method
    fix_scale=true; fix_projective=true; fix_rank=true;
    % E=G.Edges.EndNodes; node1=E(1,1); node2=E(1,2);
    % select nodes with highest degree as reference
    deg=degree(G); [~,node1]=max(deg);
    N = neighbors(G,node1); [~,ii]=max(deg(N)); node2=N(ii);
    [S_direct,~,Cams,Funds]=solvability_matrix_skew(G,fix_scale,fix_projective,fix_rank,node1,node2);
    [issolvableD,lambdaD]=direct_finite_solvability(S_direct,method);
    disp(['Direct method: Is the graph finite-solvable? ' num2str(issolvableD)])
    disp(' ')

end

%% 




