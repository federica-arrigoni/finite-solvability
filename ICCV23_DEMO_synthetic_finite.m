%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to test finite solvability on small synthetic cases

clc, clear, close all
addpath(genpath('./'))

%% Select one example

% % Graph with 9 nodes that is finite solvable but not solvable (2 solutions)
% load('./GraphsSynthetic/G9.mat'); n=size(A,1);

% % Minimal solvable graph with 90 nodes
% load('./GraphsSynthetic/G90.mat'); n=size(A,1); 

% Random graph
n=15; % nodes
A=rand(n)<0.25;
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
    lambda_Trager
    disp(['Trager et al.: Is the graph finite-solvable? ' num2str(issolvable_Trager)])

    % Our method 
    [S_reduced,nE]=solvability_matrix_simplified(G,C,true);
    [issolvable,lambda]=finite_solvability(S_reduced,nE,method);
    lambda
    disp(['Our method: Is the graph finite-solvable? ' num2str(issolvable)])

end


