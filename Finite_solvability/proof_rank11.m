
%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%% Proof of Lemma 1
 
syms c0 c1 c2 c3 real

%% Equations from Trager et al. ECCV 2018

L=sym(zeros(20,16));
L(1,14)=c2; L(1,10)=-c3;
L(2,13)=c2; L(2,9)=-c3;
L(3,15)=c1; L(3,7)=-c3;
L(4,13)=c1; L(4,5)=-c3;
L(5,12)=c1; L(5,8)=-c2;
L(6,9)=c1; L(6,5)=-c2;
L(7,15)=c0; L(7,3)=-c3;
L(8,14)=c0; L(8,2)=-c3;
L(9,12)=c0; L(9,4)=-c2;
L(10,10)=c0; L(10,2)=-c2;
L(11,8)=c0; L(11,4)=-c1;
L(12,7)=c0; L(12,3)=-c1;

L(13,11)=c1; L(13,16)=-c1; L(13,7)=-c2; L(13,8)=c3;
L(14,10)=c1; L(14,6)=-c2; L(14,16)=c2; L(14,12)=-c3;
L(15,13)=c0; L(15,15)=-c2; L(15,1)=-c3; L(15,11)=c3;
L(16,11)=c0; L(16,16)=-c0; L(16,3)=-c2; L(16,4)=c3;
L(17,9)=c0; L(17,1)=-c2; L(17,16)=c2; L(17,12)=-c3;
L(18,14)=c1; L(18,15)=-c2; L(18,6)=-c3; L(18,11)=c3;
L(19,6)=c0; L(19,16)=-c0; L(19,2)=-c1; L(19,4)=c3;
L(20,5)=c0; L(20,1)=-c1; L(20,16)=c1; L(20,8)=-c3;

size(L)
rank(L)

%% 11 independent equations (our formulation - ICCV 2023)

K=L([1:5 7:9 13:15],:);
rank(K)

% other set of 11 independent equations
U=L([1:5 7:9 15 16 18],:);
rank(U)






