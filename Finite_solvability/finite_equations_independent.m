%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
 
function A=finite_equations_independent(c) 
% INPUT: camera centre
% output: 11x16 matrix associated with the camera centre, which is the
% building block of the REDUCED finite solvability equations


% Legacy Federica's version
% A=finite_equations(c);
% % select only independent equations (11 out of 20)
% A=A([1:5 7:9 13:15],:);


C1=c(1);
C2=c(2);
C3=c(3);
C4=c(4);


i = [11
     7
     6
     8
     4
    10
     3
     9
     5
     9
     2
     1
    10
     9
    11
     5
     8
    10
     2
     4
    11
     1
     7
     3
     6
    11
     9
    10];


j = [ 1
     2
     3
     4
     5
     6
     7
     7
     8
     8
     9
    10
    10
    11
    11
    12
    12
    12
    13
    13
    13
    14
    14
    15
    15
    15
    16
    16];

 
v =[-C4
-C4
-C4
-C3
-C4
-C3
-C4
-C3
-C3
 C4
-C4
-C4
 C2
 C2
 C4
 C2
 C1
-C4
 C3
 C2
 C1
 C3
 C1
 C2
 C1
-C3
-C2
 C3];

A = sparse(i,j,v,11,16);

end