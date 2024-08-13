

function [A,m,Cams,Funds]=solvability_matrix_skew(G,fix_scale,fix_projective,fix_rank,node1,node2)
% Input:
% G = viewing graph
%
% Output:
% A = finite solvability matrix
% m = number of edges in G
% Cams = 3x4xn matrix containing cameras
% Funds = 3x3xm matrix containing fundamental matrices

if (nargin<2)
    node1=1;
    node2=2;
    fix_projective=true;
    fix_rank=true;
    fix_scale=true;
end

%% Properties of the graph

n=numnodes(G); % number of nodes in input graph
m=numedges(G); % number of edges in input graph
E=G.Edges.EndNodes; % edges in the graph

%% Sample a configuration of cameras and fundamental matrices
[Cams,Funds]=random_cams_funds(E,n,m);

%% Build the finite solvability matrix 

disp('Building the solvability matrix...')

nEq = 10*m; % number of equations
A = spalloc(nEq,12*n,nEq*12*2); % allocate memory (#rows, #columns, #nonzero)

for k=1:m % for all the edges ...

    % let us consider edge k in the input graph
    i=E(k,1); % left node
    j=E(k,2); % right node

    P=Cams(:,:,i); % left camera
    Q=Cams(:,:,j); % right camera
    F=Funds(:,:,k); % fundamental matrix

    % compute local derivatives
    [DQ,DP]=derivatives_PQ(Q,P,F); % 10 equations

    % put local derivatives in the global matrix
    A(10*k-9:10*k,[12*i-11:12*i 12*j-11:12*j])=[DP DQ]; 

end


%% fix projective ambiguity

if fix_projective
    B=spalloc(16,12*n,16); % allocate memory (#rows, #columns, #nonzero)

    % node1 is fixed (12 parameters)
    B(1:12,12*node1-11:12*node1)=speye(12);

    % 4 parameters of node 2 are fixed (16 parameters are fixed in total)
    %B(13:16,12*node2-11:12*node2)=rand(4,12); % OK 4 random linear equations in all the entries
    %B(13:16,12*node2-11:12*node2-8)=speye(4); % does not work (fixing 1st column+ 1 entry in the 2nd column)
    ind=12*node2-12; ind1_row=ind+[1 4 7 10]; B(13:16,ind1_row)=speye(4); % OK (fixing 1st row)

else
    B=[];
end

%% fix scale ambiguity

if fix_scale
    % for each camera the entries sum to 1 (or something else)
    D=sparse(n,12*n,12*n);  % allocate memory (#rows, #columns, #nonzero)
    for i=1:n
        D(i,12*i-11:12*i)=1;
    end
    D(node1,:)=[];
else
    D=[];
end

%% Include full-rank constraints

if fix_rank
    E=sparse(n,13*n,13*n);  % allocate memory (#rows, #columns, #nonzero)

    for i=1:n
        [DP,Dz]=derivatives_P_rank(Cams(:,:,i),1);
        E(i,[12*i-11:12*i 12*n+i])=[DP Dz];
    end

    E(node1,:)=[];
    E(:,12*n+node1)=[];

    A=[A sparse(size(A,1),n-1)];
    B=[B sparse(size(B,1),n-1)];
    D=[D sparse(size(D,1),n-1)];
else
    E=[];
end


%% Append all matrices

A=[A;B;D;E];


end




%%

function [DP,DQ]=derivatives_PQ(P,Q,F)

%% Diagonal Equations
%% 1st Equation

DDQ1=[F(:,1)'*P(:,1); % Derivatives of 1st equation wrt Q
    F(:,2)'*P(:,1);
    F(:,3)'*P(:,1);
    zeros(9,1)];

DDP1=[F(1,:)*Q(:,1);  % Deravitives of 1st equation wrt P
    F(2,:)*Q(:,1);
    F(3,:)*Q(:,1);
    zeros(9,1)];
% Same form as Q but we use columns instead of rows of F

%% 2nd Equation

DDQ2=[zeros(3,1);
    F(:,1)'*P(:,2); % Derivatives of 2nd equation wrt Q
    F(:,2)'*P(:,2);
    F(:,3)'*P(:,2);
    zeros(6,1)];

DDP2=[zeros(3,1);
    F(1,:)*Q(:,2);  % Derivatives of 2nd equation wrt P
    F(2,:)*Q(:,2);
    F(3,:)*Q(:,2);
    zeros(6,1)];

%% 3rd Equation

DDQ3=[zeros(6,1);
    F(:,1)'*P(:,3); % Derivatives of 3nd equation wrt Q
    F(:,2)'*P(:,3);
    F(:,3)'*P(:,3);
    zeros(3,1)];

DDP3=[zeros(6,1);
    F(1,:)*Q(:,3);  % Derivatives of 3rd equation wrt P
    F(2,:)*Q(:,3);
    F(3,:)*Q(:,3);
    zeros(3,1)];

%% 4th Equation

DDQ4=[zeros(9,1);
    F(:,1)'*P(:,4); % Derivatives of 4th equation wrt Q
    F(:,2)'*P(:,4);
    F(:,3)'*P(:,4)];

DDP4=[zeros(9,1);
    F(1,:)*Q(:,4);  % Derivatives of 4th equation wrt P
    F(2,:)*Q(:,4);
    F(3,:)*Q(:,4)];

%% Other Equations
%% 5th Equation

DDQ5=[F(:,1)'*P(:,2); % Derivatives of 5th equation wrt Q
    F(:,2)'*P(:,2);
    F(:,3)'*P(:,2);
    F(:,1)'*P(:,1);
    F(:,2)'*P(:,1);
    F(:,3)'*P(:,1);
    zeros(6,1)];

DDP5=[F(1,:)*Q(:,2); % Derivatives of 5th equation wrt P
    F(2,:)*Q(:,2);
    F(3,:)*Q(:,2);
    F(1,:)*Q(:,1);
    F(2,:)*Q(:,1);
    F(3,:)*Q(:,1);
    zeros(6,1)];

%% 6th Equation

DDQ6=[F(:,1)'*P(:,3); % Derivatives of 6th equation wrt Q
    F(:,2)'*P(:,3);
    F(:,3)'*P(:,3);
    zeros(3,1);
    F(:,1)'*P(:,1);
    F(:,2)'*P(:,1);
    F(:,3)'*P(:,1);
    zeros(3,1)];

DDP6=[F(1,:)*Q(:,3); % Derivatives of 6th equation wrt P
    F(2,:)*Q(:,3);
    F(3,:)*Q(:,3);
    zeros(3,1);
    F(1,:)*Q(:,1);
    F(2,:)*Q(:,1);
    F(3,:)*Q(:,1);
    zeros(3,1)];

%% 7th Equation

DDQ7=[F(:,1)'*P(:,4); % Derivatives of 7th equation wrt Q
    F(:,2)'*P(:,4);
    F(:,3)'*P(:,4);
    zeros(6,1);
    F(:,1)'*P(:,1);
    F(:,2)'*P(:,1);
    F(:,3)'*P(:,1)];

DDP7=[F(1,:)*Q(:,4); % Derivatives of 7th equation wrt P
    F(2,:)*Q(:,4);
    F(3,:)*Q(:,4);
    zeros(6,1);
    F(1,:)*Q(:,1);
    F(2,:)*Q(:,1);
    F(3,:)*Q(:,1)];

%% 8th Equation

DDQ8=[zeros(3,1);
    F(:,1)'*P(:,3); % Derivatives of 8th equation wrt Q
    F(:,2)'*P(:,3);
    F(:,3)'*P(:,3);
    F(:,1)'*P(:,2);
    F(:,2)'*P(:,2);
    F(:,3)'*P(:,2);
    zeros(3,1);];

DDP8=[zeros(3,1);
    F(1,:)*Q(:,3); % Derivatives of 8th equation wrt P
    F(2,:)*Q(:,3);
    F(3,:)*Q(:,3);
    F(1,:)*Q(:,2);
    F(2,:)*Q(:,2);
    F(3,:)*Q(:,2);
    zeros(3,1);];

%% 9th Equation

DDQ9=[zeros(3,1);
    F(:,1)'*P(:,4); % Derivatives of 9th equation wrt Q
    F(:,2)'*P(:,4);
    F(:,3)'*P(:,4);
    zeros(3,1);
    F(:,1)'*P(:,2);
    F(:,2)'*P(:,2);
    F(:,3)'*P(:,2)];

DDP9=[zeros(3,1);
    F(1,:)*Q(:,4); % Derivatives of 9th equation wrt P
    F(2,:)*Q(:,4);
    F(3,:)*Q(:,4);
    zeros(3,1);
    F(1,:)*Q(:,2);
    F(2,:)*Q(:,2);
    F(3,:)*Q(:,2)];

%% 10th Equation

DDQ10=[zeros(6,1);
    F(:,1)'*P(:,4); % Derivatives of 10th equation wrt Q
    F(:,2)'*P(:,4);
    F(:,3)'*P(:,4);
    F(:,1)'*P(:,3);
    F(:,2)'*P(:,3);
    F(:,3)'*P(:,3)];

DDP10=[zeros(6,1);
    F(1,:)*Q(:,4); % Derivatives of 10th equation wrt P
    F(2,:)*Q(:,4);
    F(3,:)*Q(:,4);
    F(1,:)*Q(:,3);
    F(2,:)*Q(:,3);
    F(3,:)*Q(:,3)];

%% Concatenate all equations

DP=[DDP1 DDP2 DDP3 DDP4 DDP5 DDP6 DDP7 DDP8 DDP9 DDP10]';
DQ=[DDQ1 DDQ2 DDQ3 DDQ4 DDQ5 DDQ6 DDQ7 DDQ8 DDQ9 DDQ10]';

end

%%

function [DP,Dz]=derivatives_P_rank(P,z)
% imposes that camera is full-rank by using a 4x4 matrix with random row
% and imposing that its determinant is nonzero (thank to the introduction
% of an auxiliary variable) 
% z*det([P;rand(1,4)])+1=0

a=rand(1,4);

Dz=P(1,1)*P(2,2)*P(3,3)*a(4) - P(1,1)*P(2,2)*P(3,4)*a(3) - P(1,1)*P(2,3)*P(3,2)*a(4) + P(1,1)*P(2,3)*P(3,4)*a(2) + P(1,1)*P(2,4)*P(3,2)*a(3) - P(1,1)*P(2,4)*P(3,3)*a(2) - P(1,2)*P(2,1)*P(3,3)*a(4) + P(1,2)*P(2,1)*P(3,4)*a(3) + P(1,2)*P(2,3)*P(3,1)*a(4) - P(1,2)*P(2,3)*P(3,4)*a(1) - P(1,2)*P(2,4)*P(3,1)*a(3) + P(1,2)*P(2,4)*P(3,3)*a(1) + P(1,3)*P(2,1)*P(3,2)*a(4) - P(1,3)*P(2,1)*P(3,4)*a(2) - P(1,3)*P(2,2)*P(3,1)*a(4) + P(1,3)*P(2,2)*P(3,4)*a(1) + P(1,3)*P(2,4)*P(3,1)*a(2) - P(1,3)*P(2,4)*P(3,2)*a(1) - P(1,4)*P(2,1)*P(3,2)*a(3) + P(1,4)*P(2,1)*P(3,3)*a(2) + P(1,4)*P(2,2)*P(3,1)*a(3) - P(1,4)*P(2,2)*P(3,3)*a(1) - P(1,4)*P(2,3)*P(3,1)*a(2) + P(1,4)*P(2,3)*P(3,2)*a(1);

DP=[ z*(P(2,2)*P(3,3)*a(4) - P(2,2)*P(3,4)*a(3) - P(2,3)*P(3,2)*a(4) + P(2,3)*P(3,4)*a(2) + P(2,4)*P(3,2)*a(3) - P(2,4)*P(3,3)*a(2)),...
-z*(P(1,2)*P(3,3)*a(4) - P(1,2)*P(3,4)*a(3) - P(1,3)*P(3,2)*a(4) + P(1,3)*P(3,4)*a(2) + P(1,4)*P(3,2)*a(3) - P(1,4)*P(3,3)*a(2)),...
 z*(P(1,2)*P(2,3)*a(4) - P(1,2)*P(2,4)*a(3) - P(1,3)*P(2,2)*a(4) + P(1,3)*P(2,4)*a(2) + P(1,4)*P(2,2)*a(3) - P(1,4)*P(2,3)*a(2)),...
-z*(P(2,1)*P(3,3)*a(4) - P(2,1)*P(3,4)*a(3) - P(2,3)*P(3,1)*a(4) + P(2,3)*P(3,4)*a(1) + P(2,4)*P(3,1)*a(3) - P(2,4)*P(3,3)*a(1)),...
 z*(P(1,1)*P(3,3)*a(4) - P(1,1)*P(3,4)*a(3) - P(1,3)*P(3,1)*a(4) + P(1,3)*P(3,4)*a(1) + P(1,4)*P(3,1)*a(3) - P(1,4)*P(3,3)*a(1)),...
-z*(P(1,1)*P(2,3)*a(4) - P(1,1)*P(2,4)*a(3) - P(1,3)*P(2,1)*a(4) + P(1,3)*P(2,4)*a(1) + P(1,4)*P(2,1)*a(3) - P(1,4)*P(2,3)*a(1)),...
 z*(P(2,1)*P(3,2)*a(4) - P(2,1)*P(3,4)*a(2) - P(2,2)*P(3,1)*a(4) + P(2,2)*P(3,4)*a(1) + P(2,4)*P(3,1)*a(2) - P(2,4)*P(3,2)*a(1)),...
-z*(P(1,1)*P(3,2)*a(4) - P(1,1)*P(3,4)*a(2) - P(1,2)*P(3,1)*a(4) + P(1,2)*P(3,4)*a(1) + P(1,4)*P(3,1)*a(2) - P(1,4)*P(3,2)*a(1)),...
 z*(P(1,1)*P(2,2)*a(4) - P(1,1)*P(2,4)*a(2) - P(1,2)*P(2,1)*a(4) + P(1,2)*P(2,4)*a(1) + P(1,4)*P(2,1)*a(2) - P(1,4)*P(2,2)*a(1)),...
-z*(P(2,1)*P(3,2)*a(3) - P(2,1)*P(3,3)*a(2) - P(2,2)*P(3,1)*a(3) + P(2,2)*P(3,3)*a(1) + P(2,3)*P(3,1)*a(2) - P(2,3)*P(3,2)*a(1)),...
 z*(P(1,1)*P(3,2)*a(3) - P(1,1)*P(3,3)*a(2) - P(1,2)*P(3,1)*a(3) + P(1,2)*P(3,3)*a(1) + P(1,3)*P(3,1)*a(2) - P(1,3)*P(3,2)*a(1)),...
-z*(P(1,1)*P(2,2)*a(3) - P(1,1)*P(2,3)*a(2) - P(1,2)*P(2,1)*a(3) + P(1,2)*P(2,3)*a(1) + P(1,3)*P(2,1)*a(2) - P(1,3)*P(2,2)*a(1))];

end



