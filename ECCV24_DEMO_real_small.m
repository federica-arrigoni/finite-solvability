%% Adapted from the code of the following paper:
%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to test finite solvability on small-scale real data
%% The proposed approach is compared with Trager et al. ECCV 2018

clc, clear, close all
addpath(genpath('./'))

%% SMALL Projective SfM datasets

%dataname='Dino319'; % 36 nodes, FINITE SOLVABLE
%dataname='Dino4983'; % 36 nodes, FINITE SOLVABLE
%dataname='Corridor'; % 11 nodes, complete
%dataname='House'; % 10 nodes, complete
dataname='GustavVasa'; % 18 nodes, FINITE SOLVABLE
%dataname='FolkeFilbyter'; % 40 nodes, FINITE SOLVABLE
%dataname='ParkGate'; % 34 nodes, nearly complete, FINITE SOLVABLE
%dataname='Nijo'; % 19 nodes, complete
%dataname='DrinkingFountain'; % 14 nodes, complete
%dataname='GoldenStatue'; % 18 nodes, complete
%dataname='JonasAhls'; % 40 nodes, FINITE SOLVABLE
%dataname='DeGuerre'; % 35 nodes, complete

%% MEDIUM Projective SfM datasets

%dataname='Dome'; % 85 nodes, complete
%dataname='AlcatrazCourtyard'; % n=133, nearly complete
%dataname='AlcatrazWaterTower'; % n=172, complete
%dataname='Cherub'; %avgDeg=41, n=65, FINITE SOLVABLE 
%dataname='Pumpkin'; % avgDeg=126, n=195, FINITE SOLVABLE
%dataname='Sphinx'; %avgDeg=38, n=70, FINITE SOLVABLE 
%dataname='TorontoUniversity'; %avgDeg=25, n=77, FINITE SOLVABLE
%dataname='SriThendayuthapani'; % complete, n=98
%dataname='PortaSanDonato'; %complete, n=141
%dataname='BuddahTooth'; % avgDeg=118, n=162, FINITE SOLVABLE
%dataname='TsarNikolaiI'; % avgDeg=50, n=98, FINITE SOLVABLE 
%dataname='SmolnyCathedral'; %complete, n=131
%dataname='SkansenKronan'; % avgDeg=114, n=131, nearly complete

load(['./GraphsReal/PSFM/' dataname '.mat']);

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

%% check connected and rigidity

% ncomp=graphconncomp(A);
% rigidityDec = ParallelRigidityTest(A,3);
% 
% fprintf('\n#connected components = %d \n',ncomp)
% fprintf('Is the graph rigid? %d \n\n',rigidityDec)
% if ~rigidityDec
%     rig_id = LargestMaxParRigComp(A,3);
%     A=A(rig_id,:);
%     A=A(:,rig_id);
% end

%% Check finite solvability

method='eigs'; % Default for real data
C=rand(4,ncams)*100; % random camera centres

%% Trager et al. ECCV 2018

tic
[S_Trager,nL,mL]=solvability_matrix(G,C);
tbuild_trager=toc;

tic
[issolvable_Trager,lambda_Trager]=finite_solvability(S_Trager,nL,method);
ttest_trager=toc;

assert((size(S_Trager,1)+16+m-1)==eq_Trager)
assert(size(S_Trager,2)==unk_fs)
lambda_Trager

disp(['Trager et al.: Is the graph finite-solvable? ' num2str(issolvable_Trager)])
disp(['Trager et al.: time for building finite solvability matrix ' num2str(tbuild_trager)])
disp(['Trager et al.: time for testing finite solvability ' num2str(ttest_trager)])
disp(['Trager et al.: TOTAL TIME ' num2str(ttest_trager+tbuild_trager)])
fprintf('---------------------------------------------------------\n\n')


%% Arrigoni et al. ICCV 2023

tic
[S_reduced,nE]=solvability_matrix_simplified(G,C,true);
tbuild=toc;

tic
[issolvable,lambda]=finite_solvability(S_reduced,nE,method);
ttest=toc;

assert((size(S_reduced,1)+16+m-1)==eq_ICCV23)
assert(size(S_reduced,2)==unk_fs)
lambda

disp(['ICCV23: Is the graph finite-solvable? ' num2str(issolvable)])
disp(['ICCV23: time for building finite solvability matrix ' num2str(tbuild)])
disp(['ICCV23: time for testing finite solvability ' num2str(ttest)])
disp(['ICCV23: TOTAL TIME ' num2str(ttest+tbuild)])
fprintf('---------------------------------------------------------\n\n')


%% Direct method

%method='eigs'; 
method='eigs_sm'; % Default for real data

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

disp(['Direct method: Is the graph finite-solvable? ' num2str(issolvableD)])
disp(['Direct method: time for building finite solvability matrix ' num2str(tbuildD)])
disp(['Direct method: time for testing finite solvability ' num2str(ttestD)])
disp(['Direct method: TOTAL TIME ' num2str(ttestD+tbuildD)])
fprintf('---------------------------------------------------------\n\n')


disp(['SPEED UP (time): ' num2str((ttest+tbuild)/(ttestD+tbuildD)) ' x'])
disp(['SPEED UP (nnz): ' num2str(nnz(S_reduced)/nnz(S_direct)) ' x'])

% save(['results_' dataname])


