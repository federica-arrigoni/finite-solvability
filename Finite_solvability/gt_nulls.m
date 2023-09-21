function [N,flag,bound] = gt_nulls(A,k,tol)
%NULLS  Null space of a (possibly sparse) matrix.
%
%  [N,flag,bound] = nulls(A,k,tol) attempts to compute the null space  of A.
%
%  k      An initial estimate as to the dimension of the null space.
%         Defaults to 1 if ommitted. It is only a hint, which if
%         accurate, can speed up the computation (but not dramatically).
%
%  tol    A tolerance: a vector v is considered to be a null
%         vector if norm(v) < tol.
%         Defaults to 1e-14 if ommitted.
%
%  N      A matrix whose columns are orthonormal and are null vectors of A.
%
%  flag   0 if the null space was computed successfully.
%         1 if the null space might be larger than the number of columns in N.
%
%  bound  If flag==0, this is the dimension of the null space of A.
%         If flag==1, this is an upper bound on the dimensio of the null space.
%
%  The algorithm uses normalized inverse U iteration with condition
%  estimation and normalized inverse L1*U iteration. For futher details,
%  see Gotsman and Toledo, "On the computation of null spaces of sparse
%  rectangular matrices (manuscript, 2005).
%
%  The main computational cost in this algorithm is the sparse LU
%  factorization of A. It is suitable for any matrix that Matlab can
%  factor in sparse mode.
%
%  Author: Sivan Toledo
%  Version 1.0, July 2005
%  Refactrored by Andrea Fusiello, 2023

t0 = cputime;

if ( ~issparse(A) )
    A = sparse(A);
end

[m,n] = size(A);

if m<n
    error('nulls: A must be slim')
    % added by Andrea Fusiello
end

warning off MATLAB:nearlySingularMatrix

if (nargin < 1)
    error('Not enough input arguments.');
end

if (nargin < 2)
    k = 1;
end

if (nargin < 3)
    tol = max(m,n)*eps(normest(A));
end

t1 = cputime;
[L,U,~,colP] = lu(A,'vector');

%[L,U,P,colP] = lu(A,1.0);	% changed by Andrea Fusiello

% for j=1:n
%
%     i_colP(find(colP(:,j))) = j; % inverse column permutation
% end

% A.F equivalent to above but faster
i_colP(colP) = 1:length(colP); % inverse column permutation

t2 = cputime;
fprintf(1,'Factorization time in nulls is %.2f seconds\n',t2-t1);
fprintf(1,'nulls: nnz(L) = %.01e, nnz(U) = %.01e\n',nnz(L),nnz(U));

%    for i=1:n
%     if (U(i,i)==0)
%       % U(i,i)=eps*eps*normA;
%       U(i,i)=1e-100;
%     end
%    end

% A.F equivalent to above but faster

foo = diag(U);
U = U - diag(foo);
foo(foo==0)=1e-100;
U = U + diag(foo);

if (m > n)
    L1 = L(1:n,:);
else
    L1 = L;
end

%%%
%%% first, iterate to find the null space of U
%%%

kin  = k;
kout = 0;
while ( (kout < n) )
    % disp('U iteration')
    Nu = normalizedInverseTriangularIteration(U,kin,n,i_colP);
    [Nu,kout] = filterPrefixColumns(A,Nu,tol);
    if (kout < kin)
        break;
    else
        kin = kin*2;
    end
end

dim_null_U = kout;

%%%
%%% now check whether L1 is ill conditioned
%%%

% disp('L1 iteration (for condition estimation)')
Nl = normalizedInverseTriangularIteration(L1,1,n,i_colP);
[~,kout] = filterPrefixColumns(L1,Nl,tol); % changed tolerance

if (kout == 0)
    
    %%%
    %%% L1 is not ill conditioned, so we return the null space of U
    %%%
    
    N  = Nu;
    flag  = 0;
    bound = dim_null_U;
    
    fprintf(1,'nulls: case 1, null(A)==null(U), dim(null(A))==%d\n',bound);
    
else
    
    %%%
    %%% L1 is ill conditioned, so we run inverse iteration on L1*U
    %%%
    
    fprintf(1,'nulls: L1 is ill conditioned\n');
    
    kin  = dim_null_U + 1;
    kout = 0;
    while ( (kout < n) )
        % disp('LU iteration')
        Nlu = normalizedInverseLUIteration(L1,U,kin,n,i_colP);
        [Nlu,kout] = filterPrefixColumns(L1,U*Nlu,tol);
        if (kout < kin)
            break;
        else
            kin = kin*2;
        end
    end
    
    bound = kout;
    
    %%%
    %%% There are two special cases were we can still determine the null space
    %%%
    
    if (bound == dim_null_U)
        % We got lucky
        % L1 was ill conditioned, but the null space of A is that of U
        N     = Nu;
        flag  = 0;
        fprintf(1,'nulls: case 2, null(A)==null(U), dim(null(A))==%d\n',bound);
    elseif (bound == 1 && norm( A*Nlu(:,1) ) > tol)
        % if L1*U has one null vector, and it's not a null vector
        % of A, then A has no null vectors
        N     = zeros(m,0);
        flag  = 0;
        bound = 0;
        fprintf(1,'nulls: case 3, null(A)==null(L1*U), dim(null(A))==0\n');
    else
        [N,~] = qr([ Nu Nlu ],0);
        [N,kout] = filterPrefixColumns(A,N,tol);
        
        if ( kout == bound )
            flag = 0;
            fprintf(1,'nulls: case 4, found null(A), dim(null(A))==%d\n',bound);
        else
            flag = 1;
            fprintf(1,'nulls: case 5, %d <= dim(null(A)) <= %d\n',kout,bound);
        end
        
    end
end

t3 = cputime;
fprintf(1,'Total time in nulls is %.2f seconds\n',t3-t0);

warning on MATLAB:nearlySingularMatrix

end % of main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized inverse LU iteration                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = normalizedInverseLUIteration(L,U,k,n,i_colP)

W = rand(n,k);
[W,~] = qr(W,0);

for i=1:3
    % iter = i
    % now we have an approximation W for the right singular vectors
    Tt = U' \ W;
    Q = L' \ Tt;
    [Q,~] = qr(Q,0);
    
    % we start with an approximation Q for the left singular vectors
    Y = L\Q;
    X = U \ Y;
    [W,~] = qr(X,0);
end

W = sparse(W(i_colP,:));

if (norm( isnan(W) + isinf(W) , 1 ) > 0)
    disp('nulls: ****************************************')
    disp('nulls: *** infs or nans in iverse iteration ***')
    disp('nulls: ****************************************')
    disp('nulls: *** infs or nans in iverse iteration ***')
end

end % of function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized inverse iteration                       %
% for a triangular matrix Z                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = normalizedInverseTriangularIteration(Z,k,n,i_colP)

W =rand(n,k);
[W,~] = qr(W,0);

for i=1:3
    % now we have an approximation W for the right singular vectors
    Q = Z' \ W;
    [Q,~] = qr(Q,0);
    
    % we start with an approximation Q for the left singular vectors
    X = Z \ Q;
    [W,~] = qr(X,0);
end

W = sparse(W(i_colP,:));

if (norm( isnan(W) + isinf(W) , 1 ) > 0)
    disp('nulls: ****************************************')
    disp('nulls: *** infs or nans in iverse iteration ***')
    disp('nulls: ****************************************')
end

end % of function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter only the prefix columns of W that satisfy   %
% norm( Z*W(:,k) ) <= thresh                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Wout,kout] = filterPrefixColumns(Z,W,thresh)
[~, k] = size(W);
kout=0;
while ((kout < k) && norm(Z*W(:,kout+1)) < thresh)
    kout = kout+1;
end
Wout = W(:,1:kout);

end % of function


function [L,U,Q] = luq(A,do_pivot,tol)
%  PURPOSE: calculates the following decomposition
%
%       A = L |Ubar  0 | Q
%             |0     0 |
%
%       where Ubar is a square invertible matrix
%       and matrices L, Q are invertible.
%
% ---------------------------------------------------
%  USAGE: [L,U,Q] = luq(A,do_pivot,tol)
%  INPUT:
%         A             a sparse matrix
%         do_pivot      = 1 with column pivoting
%                       = 0 without column pivoting
%         tol           uses the tolerance tol in separating zero and
%                       nonzero values
%
%   OUTPUT:
%         L,U,Q          matrices
%
%   COMMENTS:
%         based on lu decomposition
%
% Copyright  (c) Pawel Kowal (2006)
% All rights reserved
% LREM_SOLVE toolbox is available free for noncommercial academic use only.
% pkowal3@sgh.waw.pl
[n,m]                   = size(A);
if ~issparse(A)
    A                   = sparse(A);
end
%--------------------------------------------------------------------------
%       SPECIAL CASES
%--------------------------------------------------------------------------
if size(A,1)==0
    L                   = speye(n);
    U                   = A;
    Q                   = speye(m);
    return;
end
if size(A,2)==0
    L                   = speye(n);
    U                   = A;
    Q                   = speye(m);
    return;
end
%--------------------------------------------------------------------------
%       LU DECOMPOSITION
%--------------------------------------------------------------------------
if do_pivot
    [L,U,P,Q]           = lu(A);
    Q                   = Q';
else
    [L,U,P]             = lu(A);
    Q                   = speye(m);
end
p                       = size(A,1)-size(L,2);
%LL                      = [sparse(n-p,p);speye(p)];
L                       = [P'*L P(n-p+1:n,:)'];
U                       = [U;sparse(p,m)];
%--------------------------------------------------------------------------
%       FINDS ROWS WITH ZERO AND NONZERO ELEMENTS ON THE DIAGONAL
%--------------------------------------------------------------------------
if size(U,1)==1 || size(U,2)==1
    S                   = U(1,1);
else
    S                   = diag(U);
end
I                       = find(abs(S)>tol);
Jl                      = (1:n)';
Jl(I)                   = [];
Jq                      = (1:m)';
Jq(I)                   = [];
Ubar1                   = U(I,I);
Ubar2                   = U(Jl,Jq);
Qbar1                   = Q(I,:);
Lbar1                   = L(:,I);
%--------------------------------------------------------------------------
%       ELININATES NONZEZO ELEMENTS BELOW AND ON THE RIGHT OF THE
%       INVERTIBLE BLOCK OF THE MATRIX U
%
%       UPDATES MATRICES L, Q
%--------------------------------------------------------------------------
if ~isempty(I)
    Utmp                = U(I,Jq);
    X                   = Ubar1'\U(Jl,I)';
    Ubar2               = Ubar2-X'*Utmp;
    Lbar1               = Lbar1+L(:,Jl)*X';
    X                   = Ubar1\Utmp;
    Qbar1               = Qbar1+X*Q(Jq,:);
    Utmp                = [];
    X                   = [];
end
%--------------------------------------------------------------------------
%       FINDS ROWS AND COLUMNS WITH ONLY ZERO ELEMENTS
%--------------------------------------------------------------------------
I2                      = find(max(abs(Ubar2),[],2)>tol);
I5                      = find(max(abs(Ubar2),[],1)>tol);
I3                      = Jl(I2);
I4                      = Jq(I5);
Jq(I5)                  = [];
Jl(I2)                  = [];
U                       = [];
%--------------------------------------------------------------------------
%       FINDS A PART OF THE MATRIX U WHICH IS NOT IN THE REQIRED FORM
%--------------------------------------------------------------------------
A                       = Ubar2(I2,I5);
%--------------------------------------------------------------------------
%       PERFORMS LUQ DECOMPOSITION OF THE MATRIX A
%--------------------------------------------------------------------------
[L1,U1,Q1]              = luq(A,do_pivot,tol);
%--------------------------------------------------------------------------
%       UPDATES MATRICES L, U, Q
%--------------------------------------------------------------------------
Lbar2                   = L(:,I3)*L1;
Qbar2                   = Q1*Q(I4,:);
L                       = [Lbar1 Lbar2 L(:,Jl)];
Q                       = [Qbar1; Qbar2; Q(Jq,:)];
n1                      = length(I);
n2                      = length(I3);
m2                      = length(I4);
U                       = [Ubar1 sparse(n1,m-n1);sparse(n2,n1) U1 sparse(n2,m-n1-m2);sparse(n-n1-n2,m)];
end

