

clear,clc

P = sym('P',[3 4],'real'); % camera
a = sym('a',[1 4],'real'); % random values
z = sym('z','real'); % auxiliary variable

%% Determinant of 4x4 matrix with camera and random values in the 4th row

Eq=z*det([P;a])+1;
for i=1:3
    for j=1:4
        DP(i,j)=diff(Eq,P(i,j));  
    end
end
DP(:)

Dz=diff(Eq,z);  
Dz

%%
% 
% Dz=P(1,1)*P(2,2)*P(3,3)*a(4) - P(1,1)*P(2,2)*P(3,4)*a(3) - P(1,1)*P(2,3)*P(3,2)*a(4) + P(1,1)*P(2,3)*P(3,4)*a(2) + P(1,1)*P(2,4)*P(3,2)*a(3) - P(1,1)*P(2,4)*P(3,3)*a(2) - P(1,2)*P(2,1)*P(3,3)*a(4) + P(1,2)*P(2,1)*P(3,4)*a(3) + P(1,2)*P(2,3)*P(3,1)*a(4) - P(1,2)*P(2,3)*P(3,4)*a(1) - P(1,2)*P(2,4)*P(3,1)*a(3) + P(1,2)*P(2,4)*P(3,3)*a(1) + P(1,3)*P(2,1)*P(3,2)*a(4) - P(1,3)*P(2,1)*P(3,4)*a(2) - P(1,3)*P(2,2)*P(3,1)*a(4) + P(1,3)*P(2,2)*P(3,4)*a(1) + P(1,3)*P(2,4)*P(3,1)*a(2) - P(1,3)*P(2,4)*P(3,2)*a(1) - P(1,4)*P(2,1)*P(3,2)*a(3) + P(1,4)*P(2,1)*P(3,3)*a(2) + P(1,4)*P(2,2)*P(3,1)*a(3) - P(1,4)*P(2,2)*P(3,3)*a(1) - P(1,4)*P(2,3)*P(3,1)*a(2) + P(1,4)*P(2,3)*P(3,2)*a(1);
% 
% DP=[ z*(P(2,2)*P(3,3)*a(4) - P(2,2)*P(3,4)*a(3) - P(2,3)*P(3,2)*a(4) + P(2,3)*P(3,4)*a(2) + P(2,4)*P(3,2)*a(3) - P(2,4)*P(3,3)*a(2)),...
% -z*(P(1,2)*P(3,3)*a(4) - P(1,2)*P(3,4)*a(3) - P(1,3)*P(3,2)*a(4) + P(1,3)*P(3,4)*a(2) + P(1,4)*P(3,2)*a(3) - P(1,4)*P(3,3)*a(2)),...
%  z*(P(1,2)*P(2,3)*a(4) - P(1,2)*P(2,4)*a(3) - P(1,3)*P(2,2)*a(4) + P(1,3)*P(2,4)*a(2) + P(1,4)*P(2,2)*a(3) - P(1,4)*P(2,3)*a(2)),...
% -z*(P(2,1)*P(3,3)*a(4) - P(2,1)*P(3,4)*a(3) - P(2,3)*P(3,1)*a(4) + P(2,3)*P(3,4)*a(1) + P(2,4)*P(3,1)*a(3) - P(2,4)*P(3,3)*a(1)),...
%  z*(P(1,1)*P(3,3)*a(4) - P(1,1)*P(3,4)*a(3) - P(1,3)*P(3,1)*a(4) + P(1,3)*P(3,4)*a(1) + P(1,4)*P(3,1)*a(3) - P(1,4)*P(3,3)*a(1)),...
% -z*(P(1,1)*P(2,3)*a(4) - P(1,1)*P(2,4)*a(3) - P(1,3)*P(2,1)*a(4) + P(1,3)*P(2,4)*a(1) + P(1,4)*P(2,1)*a(3) - P(1,4)*P(2,3)*a(1)),...
%  z*(P(2,1)*P(3,2)*a(4) - P(2,1)*P(3,4)*a(2) - P(2,2)*P(3,1)*a(4) + P(2,2)*P(3,4)*a(1) + P(2,4)*P(3,1)*a(2) - P(2,4)*P(3,2)*a(1)),...
% -z*(P(1,1)*P(3,2)*a(4) - P(1,1)*P(3,4)*a(2) - P(1,2)*P(3,1)*a(4) + P(1,2)*P(3,4)*a(1) + P(1,4)*P(3,1)*a(2) - P(1,4)*P(3,2)*a(1)),...
%  z*(P(1,1)*P(2,2)*a(4) - P(1,1)*P(2,4)*a(2) - P(1,2)*P(2,1)*a(4) + P(1,2)*P(2,4)*a(1) + P(1,4)*P(2,1)*a(2) - P(1,4)*P(2,2)*a(1)),...
% -z*(P(2,1)*P(3,2)*a(3) - P(2,1)*P(3,3)*a(2) - P(2,2)*P(3,1)*a(3) + P(2,2)*P(3,3)*a(1) + P(2,3)*P(3,1)*a(2) - P(2,3)*P(3,2)*a(1)),...
%  z*(P(1,1)*P(3,2)*a(3) - P(1,1)*P(3,3)*a(2) - P(1,2)*P(3,1)*a(3) + P(1,2)*P(3,3)*a(1) + P(1,3)*P(3,1)*a(2) - P(1,3)*P(3,2)*a(1)),...
% -z*(P(1,1)*P(2,2)*a(3) - P(1,1)*P(2,3)*a(2) - P(1,2)*P(2,1)*a(3) + P(1,2)*P(2,3)*a(1) + P(1,3)*P(2,1)*a(2) - P(1,3)*P(2,2)*a(1))];
