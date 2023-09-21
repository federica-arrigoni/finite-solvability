%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
 
function A=finite_equations(c) 
% INPUT: camera centre
% output: 20x16 matrix associated with the camera centre, which is the
% building block of the finite solvability equations


% Legacy Federica's version
% c0=c(1);
% c1=c(2);
% c2=c(3);
% c3=c(4);
% 
% A=sparse(20,16);
% A(1,14)=c2; A(1,10)=-c3;
% A(2,13)=c2; A(2,9)=-c3;
% A(3,15)=c1; A(3,7)=-c3;
% A(4,13)=c1; A(4,5)=-c3;
% A(5,12)=c1; A(5,8)=-c2;
% A(6,9)=c1; A(6,5)=-c2;
% A(7,15)=c0; A(7,3)=-c3;
% A(8,14)=c0; A(8,2)=-c3;
% A(9,12)=c0; A(9,4)=-c2;
% A(10,10)=c0; A(10,2)=-c2;
% A(11,8)=c0; A(11,4)=-c1;
% A(12,7)=c0; A(12,3)=-c1;
% 
% A(13,11)=c1; A(13,16)=-c1; A(13,7)=-c2; A(13,8)=c3;
% A(14,10)=c1; A(14,6)=-c2; A(14,16)=c2; A(14,12)=-c3;
% A(15,13)=c0; A(15,15)=-c2; A(15,1)=-c3; A(15,11)=c3;
% A(16,11)=c0; A(16,16)=-c0; A(16,3)=-c2; A(16,4)=c3;
% A(17,9)=c0; A(17,1)=-c2; A(17,16)=c2; A(17,12)=-c3;
% A(18,14)=c1; A(18,15)=-c2; A(18,6)=-c3; A(18,11)=c3;
% A(19,6)=c0; A(19,16)=-c0; A(19,2)=-c1; A(19,4)=c3;
% A(20,5)=c0; A(20,1)=-c1; A(20,16)=c1; A(20,8)=-c3;


C1=c(1);
C2=c(2);
C3=c(3);
C4=c(4);


i =[15
    17
    20
     8
    10
    19
     7
    12
    16
     9
    11
    16
    19
     4
     6
    20
    14
    18
    19
     3
    12
    13
     5
    11
    13
    20
     2
     6
    17
     1
    10
    14
    13
    15
    16
    18
     5
     9
    14
    17
     2
     4
    15
     1
     8
    18
     3
     7
    15
    18
    13
    14
    16
    17
    19
    20];

j =[ 1
     1
     1
     2
     2
     2
     3
     3
     3
     4
     4
     4
     4
     5
     5
     5
     6
     6
     6
     7
     7
     7
     8
     8
     8
     8
     9
     9
     9
    10
    10
    10
    11
    11
    11
    11
    12
    12
    12
    12
    13
    13
    13
    14
    14
    14
    15
    15
    15
    15
    16
    16
    16
    16
    16
    16];    

v = [-C4
-C3
-C2
-C4
-C3
-C2
-C4
-C2
-C3
-C3
-C2
 C4
 C4
-C4
-C3
 C1
-C3
-C4
 C1
-C4
 C1
-C3
-C3
 C1
 C4
-C4
-C4
 C2
 C1
-C4
 C1
 C2
 C2
 C4
 C1
 C4
 C2
 C1
-C4
-C4
 C3
 C2
 C1
 C3
 C1
 C2
 C2
 C1
-C3
-C3
-C2
 C3
-C1
 C3
-C1
 C2];

A = sparse(i,j,v,20,16);

end