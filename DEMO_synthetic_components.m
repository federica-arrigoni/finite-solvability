%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to compute maximal finite-solvable components on small synthetic cases

clc, clear, close all
addpath(genpath('./'))

%% Select one example

% % Toy examples with known ground-truth components
% load('./GraphsSynthetic/house.mat'); n=size(A,1);
% load('./GraphsSynthetic/triangles2.mat'); n=size(A,1);
% load('./GraphsSynthetic/triangles3.mat'); n=size(A,1);

% Random connected graph        
n=20; % 20 nodes
prob_edges=0.2; % 0.2 edges
A=rand(n)<prob_edges;
A=triu(A,1); A=A+A';
cc=conncomp(graph(A));
while max(cc)~=1
    A=rand(n)<prob_edges;
    A=triu(A,1); A=A+A';
    cc=conncomp(graph(A));
end

G=graph(A);
% figure, plot(G)

%% Compute the reduced finite solvability matrix

C=rand(4,n)*100; % random camera centres
[S_reduced,nE]=solvability_matrix_simplified(G,C,true);

%% Extract all maximal finite-solvable components

method = 'svd'; % use Matlab function null (based on svd) --> SMALL graphs
%method = 'nulls'; % use  Gotsman-Toledo nulls --> LARGE sparse graphs

[fcomp,U]=finite_solvable_components(S_reduced,nE,method);

% Print edges, ground-truth components (if available), found components
if exist('gt_comp','var')
    [G.Edges{:,1} gt_comp fcomp] % edges - ground-truth components - found components
else
    [G.Edges{:,1} fcomp] % edges - ground-truth components - found components
end

%% Plot all components: each (random) colour corresponds to one component

ncomp=max(fcomp); % number of components
col=rand(ncomp,3); % colors to use
col_edges=nan(nE,3); % colors assigned to edges
for k=1:nE
    col_edges(k,:)=col(fcomp(k),:);
end
G.Edges.EdgeColors=col_edges;

figure, plot(G,'EdgeColor',G.Edges.EdgeColors,'LineWidth',5,'MarkerSize',10,'NodeFontSize',15,'NodeColor','k')
title(['#components = ' num2str(ncomp)])
set(gca,'FontSize',20)

% %% Plot the largest component
% 
% largest=mode(fcomp); % index of largest component
% edges_largest=find(fcomp==largest); % edges in the largest component
% 
% E=G.Edges{:,1}(edges_largest,:);
% Glargest=graph();
% Glargest=addedge(Glargest,E(:,1),E(:,2));
% 
% figure, plot(Glargest,'EdgeColor',col(fcomp(largest),:),'LineWidth',5,'MarkerSize',10,'NodeFontSize',15,'NodeColor','k')
% title('Largest component')
% set(gca,'FontSize',20)

