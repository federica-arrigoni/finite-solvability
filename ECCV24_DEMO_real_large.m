%% Adapted from the code of the following paper:
%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to compute maximal finite-solvable components on large-scale real data

clc, clear, close all
addpath(genpath('./'))

%% 1DSfM datasets

%dataname = 'Alamo'; 
%dataname = 'Ellis_Island'; 
%dataname = 'Madrid_Metropolis'; 
%dataname = 'Montreal_Notre_Dame'; 
%dataname = 'Notre_Dame'; 
%dataname = 'NYC_Library'; 
%dataname = 'Piazza_del_Popolo'; 
%dataname = 'Piccadilly'; 
%dataname = 'Roman_Forum'; 
dataname = 'Tower_of_London'; 
%dataname = 'Union_Square'; 
%dataname = 'Vienna_Cathedral'; 
%dataname = 'Yorkminster';
%dataname = 'Gendarmenmarkt'; 
%dataname = 'Quad'; 
%dataname = 'Trafalgar'; % crashed

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

%% Compute the size of the solvability matrices

% The number of equations comprise also equations for fixing ambiguity,
% which are not counted in ICCV23 paper

nL=m; % nodes in the line graph = edges in input graph
mL=sum(deg.^2)/2-m; % edges in the line graph

unk_fs=16*nL;
eq_Trager=20*mL+16+nL-1;
eq_ICCV23=sum(11*(degree(G)-1))+16+nL-1;

unk_our=12*ncams+ncams-1;
eq_our=10*m+16+ncams-1+ncams-1;

fprintf('#unknowns ICCV23 and Trager et al.= %d \n',unk_fs)
fprintf('#unknowns our method= %d \n',unk_our)

fprintf('#equations Trager et al. = %d \n',eq_Trager)
fprintf('#equations ICCV23 = %d \n',eq_ICCV23)
fprintf('#equations our method = %d \n',eq_our)

fprintf('Ratio between #equations by Trager et al. and ICCV23 = %.1f \n',eq_Trager/eq_ICCV23)
fprintf('Ratio between #equations by ICCV23 and ours = %.1f \n',eq_ICCV23/eq_our)
fprintf('Ratio between #unknowns by ICCV23 and ours = %.1f \n',unk_fs/unk_our)
fprintf('---------------------------------------------------------\n\n')


%% Check finite solvability
%% ICCV 2023

method='eigs'; % Default for real data
tic
C=rand(4,ncams)*100; % random camera centres
[S_reduced,nE]=solvability_matrix_simplified(G,C,true);
tbuild=toc;

tic
[issolvable,lambda]=finite_solvability(S_reduced,nE,method);
ttest=toc;

assert((size(S_reduced,1)+16+m-1)==eq_ICCV23)
assert(size(S_reduced,2)==unk_fs)
lambda

%% ICCV 2023 - Extract all maximal finite-solvable components

if ~issolvable
    method = 'nulls'; % use  Gotsman-Toledo nulls --> LARGE sparse graphs

    tic
    [fcomp,U]=finite_solvable_components(S_reduced,nE,method);
    tcomp=toc;
    ncomp=max(fcomp); % number of components

else
    fcomp=ones(nE,1);
    ncomp=1;
    tcomp=0;
end

disp(' ')
disp(['ICCV23: Is the graph finite-solvable? ' num2str(issolvable)])
disp(['ICCV23: time for building finite solvability matrix ' num2str(tbuild)])
disp(['ICCV23: time for testing finite solvability ' num2str(ttest)])
disp(['ICCV23: TOTAL TESTING TIME ' num2str(ttest+tbuild)])
disp(['nnz = ' num2str(nnz(S_reduced))])
disp(['#components = ' num2str(ncomp)])
disp(['ICCV23: time for extracting components ' num2str(tcomp)])
disp(['ICCV23: TOTAL TIME ' num2str(ttest+tbuild+tcomp)])
fprintf('---------------------------------------------------------\n\n')

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

%% Direct method

method='eigs_sm'; % Default for real data
nE = numedges(G);

tic
fix_scale=true; fix_projective=true; fix_rank=true;
% E=G.Edges.EndNodes; node1=E(1,1); node2=E(1,2);
% select nodes with highest degree as reference
deg=degree(G); [~,node1]=max(deg); 
N = neighbors(G,node1); [~,ii]=max(deg(N)); node2=N(ii);
[S_direct,~,Cams,Funds]=solvability_matrix_skew(G,fix_scale,fix_projective,fix_rank,node1,node2);
tbuildD=toc;

tic
[issolvableD,lambdaD]=direct_finite_solvability(S_direct,method);
ttestD=toc;

lambdaD
assert(size(S_direct,1)==eq_our)
assert(size(S_direct,2)==unk_our)


%% Direct method - Extract all maximal finite-solvable components

if ~issolvableD
    method = 'nulls'; % use  Gotsman-Toledo nulls --> LARGE sparse graphs

    tic
    [fcompD,UD]=fast_direct_solvable_components(S_direct,G,nE,ncams,method,node1,node2);
    tcompD=toc;
    ncompD=max(fcompD); % number of components
    
else
    fcompD=ones(nE,1);
    ncompD=1;
    tcompD=0;
end

disp(' ')
disp(['Direct method: Is the graph finite-solvable? ' num2str(issolvableD)])
disp(['Direct method: time for building finite solvability matrix ' num2str(tbuildD)])
disp(['Direct method: time for testing finite solvability ' num2str(ttestD)])
disp(['Direct method: TOTAL TESTING TIME ' num2str(ttestD+tbuildD)])
disp(['nnz = ' num2str(nnz(S_direct))])
disp(['#components = ' num2str(ncompD)])
disp(['Direct method: time for extracting components ' num2str(tcompD)])
disp(['Direct method: TOTAL TIME ' num2str(ttestD+tbuildD+tcompD)])
fprintf('---------------------------------------------------------\n\n')


%% Comparison between the two methods in time and components

disp(['SPEED UP (time): ' num2str((ttest+tbuild+tcomp)/(ttestD+tbuildD+tcompD)) ' x'])
disp(['SPEED UP (nnz): ' num2str(nnz(S_reduced)/nnz(S_direct)) ' x'])

fcompD=bestMap(fcomp,fcompD);
missrate = sum(fcomp ~= fcompD) / length(fcompD);
disp(['Difference in components between our method and ICCV23: ' num2str(missrate*100) ' %'])

% save(['results_' dataname])


