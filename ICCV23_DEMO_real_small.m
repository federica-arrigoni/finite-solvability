%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% DEMO to test finite solvability on small-scale real data
%% The proposed approach is compared with Trager et al. ECCV 2018

clc, clear, close all
addpath(genpath('./'))

%% SMALL Projective SfM datasets

dataname='Dino319'; % 36 nodes, FINITE SOLVABLE
%dataname='Dino4983'; % 36 nodes, FINITE SOLVABLE
%dataname='Corridor'; % 11 nodes, complete
%dataname='House'; % 10 nodes, complete
%dataname='GustavVasa'; % 18 nodes, FINITE SOLVABLE
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
figure,plot(G)
figure,spy(A+eye(ncams))

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

% Our method
tic
[S_reduced,nE]=solvability_matrix_simplified(G,C,true);
tbuild=toc;
tic
[issolvable,lambda]=finite_solvability(S_reduced,nE,method);
ttest=toc;
lambda
disp(['Our method: Is the graph finite-solvable? ' num2str(issolvable)])
disp(['Our method: time for building finite solvability matrix ' num2str(tbuild)])
disp(['Our method: time for testing finite solvability ' num2str(ttest)])
fprintf('---------------------------------------------------------\n\n')

% Trager et al. ECCV 2018
tic
[S_Trager,nL,mL]=solvability_matrix(G,C);
tbuild_trager=toc;
tic
[issolvable_Trager,lambda_Trager]=finite_solvability(S_Trager,nL,method);
ttest_trager=toc;
lambda_Trager
disp(['Trager et al.: Is the graph finite-solvable? ' num2str(issolvable_Trager)])
disp(['Trager et al.: time for building finite solvability matrix ' num2str(tbuild_trager)])
disp(['Trager et al.: time for testing finite solvability ' num2str(ttest_trager)])
fprintf('---------------------------------------------------------\n\n')


% save(['results_' dataname])


