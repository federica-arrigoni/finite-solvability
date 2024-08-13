
%% Adapted from the code of the following paper:
%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to compute maximal finite-solvable components on small synthetic cases

clc, clear, %close all
addpath(genpath('./'))

%% Select one example

% % Toy examples with known ground-truth components
%load('./GraphsSynthetic/house.mat'); n=size(A,1);
%load('./GraphsSynthetic/triangles2.mat'); n=size(A,1);
%load('./GraphsSynthetic/triangles3.mat'); n=size(A,1);

%A(1,2)=1; A(2,3)=1; A(3,1)=1; n=3; A=A+A'; % triangle
%A(1,2)=1; A(2,3)=1; A(3,1)=1; A(1,4)=1; A(4,2)=1; n=4; A=A+A'; % square with chord
% n=3; A=ones(n); A=triu(A,1); A=A+A'; % complete graph 

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
%figure, plot(G)


%% ICCV 2023 Extract all maximal finite-solvable components

fprintf('---------------------------------------------------------\n\n')
disp('ICCV 2023')
% Compute the reduced finite solvability matrix
C=rand(4,n)*100; % random camera centres
[S_reduced,nE]=solvability_matrix_simplified(G,C,true);

size(S_reduced)

method = 'svd'; % use Matlab function null (based on svd) --> SMALL graphs
%method = 'nulls'; % use  Gotsman-Toledo nulls --> LARGE sparse graphs

[fcomp,U]=finite_solvable_components(S_reduced,nE,method);
ncomp=max(fcomp);

%% Plot all components: each (random) colour corresponds to one component

if ncomp<=7
    col=lines(ncomp);
else
    col=rand(ncomp,3); % colors to use
end
col_edges=nan(nE,3); % colors assigned to edges
for k=1:nE
    col_edges(k,:)=col(fcomp(k),:);
end
G.Edges.EdgeColors=col_edges;

figure, plot(G,'EdgeColor',G.Edges.EdgeColors,'LineWidth',4,'MarkerSize',12,'NodeFontSize',18,'NodeColor','k')
%title(['#components = ' num2str(ncomp)])
title([num2str(ncomp) ' components'])
set(gca,'FontSize',20)
grid on 

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


%% Direct Method: Extract all maximal finite-solvable components

fprintf('---------------------------------------------------------\n\n')
disp('Direct Method')

% Compute the solvability matrix
% select nodes with highest degree as reference
% E=G.Edges.EndNodes; node1=E(1,1); node2=E(1,2);
fix_scale=true; fix_projective=true; fix_rank=true;
deg=degree(G); [~,node1]=max(deg); 
N = neighbors(G,node1); [~,ii]=max(deg(N)); node2=N(ii);
[S_direct,~,Cams,Funds]=solvability_matrix_skew(G,fix_scale,fix_projective,fix_rank,node1,node2);

size(S_direct)

method = 'svd'; % use Matlab function null (based on svd) --> SMALL graphs
%method = 'nulls'; % use  Gotsman-Toledo nulls --> LARGE sparse graphs

% Compute components
verbose=false;
[fcompD,UD]=fast_direct_solvable_components(S_direct,G,nE,n,method,node1,node2,verbose);
ncompD=max(fcompD); % number of components

%% Analysis on nullspace

% figure,imagesc(UD),colorbar
% size(UD)
% 
% V=UD; V(abs(V)<1e-5)=0; figure,spy(V)
% % 
% i=1;U1=U(12*i-11:12*i,:);
% i=2;U2=U(12*i-11:12*i,:);
% i=3;U3=U(12*i-11:12*i,:);
% i=4;U4=U(12*i-11:12*i,:);
% i=5;U5=U(12*i-11:12*i,:);
% i=6;U6=U(12*i-11:12*i,:);
% i=7;U7=U(12*i-11:12*i,:);

%% Print statistics and results

fprintf('---------------------------------------------------------\n\n')
fcompD=bestMap(fcomp,fcompD); % align labels
missrate = sum(fcomp ~= fcompD) / length(fcompD);

% Print edges, ground-truth components (if available), found components
if exist('gt_comp','var')
    disp(['Edge (left/right node) - Ground-truth - ICCV23 - Our method'])
    [G.Edges{:,1} gt_comp fcomp fcompD] % edges - ground-truth components - found components
else
    disp(['Edge (left/right node) - ICCV23 - Our method'])
    [G.Edges{:,1} fcomp fcompD] % edges - ground-truth components - found components
end

disp(['ICCV 2023: number of components = ' num2str(ncomp)])
disp(['Direct Method: number of components = ' num2str(ncompD)])
disp(['Difference in components between our method and ICCV23: ' num2str(missrate*100) ' %'])


