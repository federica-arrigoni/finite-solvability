%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to compute maximal finite-solvable components on large-scale real data

clc, clear, close all
addpath(genpath('./'))

%% 1DSfM datasets

%dataname = 'Alamo'; 
dataname = 'Ellis_Island'; 
%dataname = 'Madrid_Metropolis'; 
%dataname = 'Montreal_Notre_Dame'; 
%dataname = 'Notre_Dame'; 
%dataname = 'NYC_Library'; 
%dataname = 'Piazza_del_Popolo'; 
%dataname = 'Piccadilly'; 
%dataname = 'Roman_Forum'; 
%dataname = 'Tower_of_London'; 
%dataname = 'Union_Square'; 
%dataname = 'Vienna_Cathedral'; 
%dataname = 'Yorkminster';
%dataname = 'Gendarmenmarkt'; 
%dataname = 'Quad';  
%dataname = 'Trafalgar'; 

load(['./GraphsReal/1DSFM_rigid/' dataname '.mat']); % these graphs are already rigid

%% GRAPH properties

ncams=size(A,1);
G=graph(A);
%figure,plot(G)
%figure,spy(A+eye(ncams))

m_minimal=ceil((11*ncams-15)/7); % edges in minimal graph
fraction_edges=100*nnz(triu(A,1))/nchoosek(ncams,2); % percentage of edges
deg=sum(A,2); % degree of nodes
m=nnz(triu(A,1)); % edges

fprintf('#nodes = %d \n',ncams)
fprintf('#edges = %d \n',m)
fprintf('There are %.2f %% edges\n\n',fraction_edges)

%% Compute the size of the solvability matrix

nL=m; % nodes in the line graph
mL=sum(deg.^2)/2-m; % edges in the line graph

unk_fs=16*nL;
fprintf('#unknowns = %d \n',unk_fs)

eq_fs=20*mL;
fprintf('#equations Trager et al. = %d \n',eq_fs)

eq_simple_independent=sum(11*(degree(G)-1));
fprintf('#equations our formulation = %d \n',eq_simple_independent)

fprintf('Ratio between #equations by Trager et al. and ours = %.1f \n',eq_fs/eq_simple_independent)
fprintf('---------------------------------------------------------\n\n')


%% Check finite solvability

method='eigs'; % Default for real data
C=rand(4,ncams)*100; % random camera centres

% Our method
tic
[S_reduced,nE]=solvability_matrix_simplified(G,C,true);
tbuild=toc;
disp(['Our method: time for building finite solvability matrix ' num2str(tbuild)])

tic
[issolvable,lambda]=finite_solvability(S_reduced,nE,method);
ttest=toc;
lambda

disp(['Our method: time for testing finite solvability ' num2str(ttest)])
disp(['Our method: Is the graph finite-solvable? ' num2str(issolvable)])

% save(['results_' dataname])

%% Extract all maximal finite-solvable components

if ~issolvable
    method = 'nulls'; % use  Gotsman-Toledo nulls --> LARGE sparse graphs

    tic
    [fcomp,U]=finite_solvable_components(S_reduced,nE,method);
    tcomp=toc;
    disp(['Our method: time for extracting components ' num2str(tcomp)])

    ncomp=max(fcomp); % number of components
    disp(['#components = ' num2str(ncomp)])
else
    fcomp=ones(nE,1);
    ncomp=1;
end

%% Plot all components: each (random) colour corresponds to one component

% col=rand(ncomp,3); % colors to use
% col_edges=nan(nE,3); % colors assigned to edges
% for k=1:nE
%     col_edges(k,:)=col(fcomp(k),:);
% end
% G.Edges.EdgeColors=col_edges;
% 
% figure, plot(G,'EdgeColor',G.Edges.EdgeColors,'LineWidth',5,'MarkerSize',10,'NodeFontSize',15,'NodeColor','k')
% title(['#components = ' num2str(ncomp)])
% set(gca,'FontSize',20)






