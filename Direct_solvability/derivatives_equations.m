

P = sym('P',[3 4],'real'); % right camera
Q = sym('Q',[3 4],'real'); % left camera
F = sym('F',[3 3],'real'); % fundamental matrix

S=P'*F*Q; % this matrix should be skew-symmetric

% diagonal elements should be zero
Eq1=S(1,1)
Eq2=S(2,2)
Eq3=S(3,3)
Eq4=S(4,4)

% off-diagonal elements should be the opposite
Eq5=S(1,2)+S(2,1)
Eq6=S(1,3)+S(3,1)
Eq7=S(1,4)+S(4,1)
Eq8=S(2,3)+S(3,2)
Eq9=S(2,4)+S(4,2)
Eq10=S(3,4)+S(4,3)

%% diagonal equations

for i=1:3
    for j=1:4
        DP1(i,j)=diff(Eq1,P(i,j));
        DP2(i,j)=diff(Eq2,P(i,j));
        DP3(i,j)=diff(Eq3,P(i,j));
        DP4(i,j)=diff(Eq4,P(i,j));

        DQ1(i,j)=diff(Eq1,Q(i,j));
        DQ2(i,j)=diff(Eq2,Q(i,j));
        DQ3(i,j)=diff(Eq3,Q(i,j));
        DQ4(i,j)=diff(Eq4,Q(i,j));        
    end
end

%% other equations

for i=1:3
    for j=1:4
        DP5(i,j)=diff(Eq5,P(i,j));
        DP6(i,j)=diff(Eq6,P(i,j));
        DP7(i,j)=diff(Eq7,P(i,j));
        DP8(i,j)=diff(Eq8,P(i,j));
        DP9(i,j)=diff(Eq9,P(i,j));
        DP10(i,j)=diff(Eq10,P(i,j));  

        DQ5(i,j)=diff(Eq5,Q(i,j));
        DQ6(i,j)=diff(Eq6,Q(i,j));
        DQ7(i,j)=diff(Eq7,Q(i,j));
        DQ8(i,j)=diff(Eq8,Q(i,j));
        DQ9(i,j)=diff(Eq9,Q(i,j));
        DQ10(i,j)=diff(Eq10,Q(i,j));  
    end
end

%% SIMPLIFIED FORMULAS for derivatives
%% Diagonal Equations
%% 1st Equation

DDQ1=[F(:,1)'*P(:,1); % Derivatives of 1st equation wrt Q
F(:,2)'*P(:,1);
F(:,3)'*P(:,1);
zeros(9,1)];

% check that simplified formulas coincide with symbolic ones
assert(all(double(DQ1(:)-DDQ1(:)))==0);

DDP1=[F(1,:)*Q(:,1);  % Deravitives of 1st equation wrt P
F(2,:)*Q(:,1);
F(3,:)*Q(:,1);
zeros(9,1)];
% Same form as Q but we use columns instead of rows of F

% check that simplified formulas coincide with symbolic ones
assert(all(double(DP1(:)-DDP1(:)))==0);

%% 2nd Equation

DDQ2=[zeros(3,1);
    F(:,1)'*P(:,2); % Derivatives of 2nd equation wrt Q
    F(:,2)'*P(:,2);
    F(:,3)'*P(:,2);
    zeros(6,1)];
assert(all(double(DQ2(:)-DDQ2(:)))==0); % check 

DDP2=[zeros(3,1);
    F(1,:)*Q(:,2);  % Derivatives of 2nd equation wrt P
    F(2,:)*Q(:,2);
    F(3,:)*Q(:,2);
    zeros(6,1)];
assert(all(double(DP2(:)-DDP2(:)))==0); % check

%% 3rd Equation

DDQ3=[zeros(6,1);
    F(:,1)'*P(:,3); % Derivatives of 3nd equation wrt Q
    F(:,2)'*P(:,3);
    F(:,3)'*P(:,3);
    zeros(3,1)];
assert(all(double(DQ3(:)-DDQ3(:)))==0); % check 

DDP3=[zeros(6,1);
    F(1,:)*Q(:,3);  % Derivatives of 3rd equation wrt P
    F(2,:)*Q(:,3);
    F(3,:)*Q(:,3);
    zeros(3,1)];
assert(all(double(DP3(:)-DDP3(:)))==0); % check

%% 4th Equation

DDQ4=[zeros(9,1);
    F(:,1)'*P(:,4); % Derivatives of 4th equation wrt Q
    F(:,2)'*P(:,4);
    F(:,3)'*P(:,4)];
assert(all(double(DQ4(:)-DDQ4(:)))==0); % check 

DDP4=[zeros(9,1);
    F(1,:)*Q(:,4);  % Derivatives of 4th equation wrt P
    F(2,:)*Q(:,4);
    F(3,:)*Q(:,4)];
assert(all(double(DP4(:)-DDP4(:)))==0); % check

%% Other Equations
%% 5th Equation

DDQ5=[F(:,1)'*P(:,2); % Derivatives of 5th equation wrt Q
F(:,2)'*P(:,2); 
F(:,3)'*P(:,2);
F(:,1)'*P(:,1);
F(:,2)'*P(:,1);
F(:,3)'*P(:,1);
zeros(6,1)];
assert(all(double(DQ5(:)-DDQ5(:)))==0); % check

DDP5=[F(1,:)*Q(:,2); % Derivatives of 5th equation wrt P
F(2,:)*Q(:,2); 
F(3,:)*Q(:,2);
F(1,:)*Q(:,1);
F(2,:)*Q(:,1);
F(3,:)*Q(:,1);
zeros(6,1)];
assert(all(double(DP5(:)-DDP5(:)))==0); % check 

%% 6th Equation

DDQ6=[F(:,1)'*P(:,3); % Derivatives of 6th equation wrt Q
F(:,2)'*P(:,3); 
F(:,3)'*P(:,3);
zeros(3,1);
F(:,1)'*P(:,1);
F(:,2)'*P(:,1);
F(:,3)'*P(:,1);
zeros(3,1)];
assert(all(double(DQ6(:)-DDQ6(:)))==0); % check

DDP6=[F(1,:)*Q(:,3); % Derivatives of 6th equation wrt P
F(2,:)*Q(:,3); 
F(3,:)*Q(:,3);
zeros(3,1);
F(1,:)*Q(:,1);
F(2,:)*Q(:,1);
F(3,:)*Q(:,1);
zeros(3,1)];
assert(all(double(DP6(:)-DDP6(:)))==0); % check 

%% 7th Equation

DDQ7=[F(:,1)'*P(:,4); % Derivatives of 7th equation wrt Q
F(:,2)'*P(:,4); 
F(:,3)'*P(:,4);
zeros(6,1);
F(:,1)'*P(:,1);
F(:,2)'*P(:,1);
F(:,3)'*P(:,1)];
assert(all(double(DQ7(:)-DDQ7(:)))==0); % check

DDP7=[F(1,:)*Q(:,4); % Derivatives of 7th equation wrt P
F(2,:)*Q(:,4); 
F(3,:)*Q(:,4);
zeros(6,1);
F(1,:)*Q(:,1);
F(2,:)*Q(:,1);
F(3,:)*Q(:,1)];
assert(all(double(DP7(:)-DDP7(:)))==0); % check 

%% 8th Equation

DDQ8=[zeros(3,1);
F(:,1)'*P(:,3); % Derivatives of 8th equation wrt Q
F(:,2)'*P(:,3); 
F(:,3)'*P(:,3);
F(:,1)'*P(:,2);
F(:,2)'*P(:,2);
F(:,3)'*P(:,2);
zeros(3,1);];
assert(all(double(DQ8(:)-DDQ8(:)))==0); % check

DDP8=[zeros(3,1);
F(1,:)*Q(:,3); % Derivatives of 8th equation wrt P
F(2,:)*Q(:,3); 
F(3,:)*Q(:,3);
F(1,:)*Q(:,2);
F(2,:)*Q(:,2);
F(3,:)*Q(:,2);
zeros(3,1);];
assert(all(double(DP8(:)-DDP8(:)))==0); % check 

%% 9th Equation

DDQ9=[zeros(3,1);
F(:,1)'*P(:,4); % Derivatives of 9th equation wrt Q
F(:,2)'*P(:,4); 
F(:,3)'*P(:,4);
zeros(3,1);
F(:,1)'*P(:,2);
F(:,2)'*P(:,2);
F(:,3)'*P(:,2)];
assert(all(double(DQ9(:)-DDQ9(:)))==0); % check

DDP9=[zeros(3,1);
F(1,:)*Q(:,4); % Derivatives of 9th equation wrt P
F(2,:)*Q(:,4); 
F(3,:)*Q(:,4);
zeros(3,1);
F(1,:)*Q(:,2);
F(2,:)*Q(:,2);
F(3,:)*Q(:,2)];
assert(all(double(DP9(:)-DDP9(:)))==0); % check 

%% 10th Equation

DDQ10=[zeros(6,1);
F(:,1)'*P(:,4); % Derivatives of 10th equation wrt Q
F(:,2)'*P(:,4); 
F(:,3)'*P(:,4);
F(:,1)'*P(:,3);
F(:,2)'*P(:,3);
F(:,3)'*P(:,3)];
assert(all(double(DQ10(:)-DDQ10(:)))==0); % check

DDP10=[zeros(6,1);
F(1,:)*Q(:,4); % Derivatives of 10th equation wrt P
F(2,:)*Q(:,4); 
F(3,:)*Q(:,4);
F(1,:)*Q(:,3);
F(2,:)*Q(:,3);
F(3,:)*Q(:,3)];
assert(all(double(DP10(:)-DDP10(:)))==0); % check 






