%% Generate random cameras and corresponding fundamental matrices

function [Cams,Funds]=random_cams_funds(E,n,m)

% generate random full-rank cameras
Cams=zeros(3,4,n);
sc=10; % define a range for the entries of cameras/fundamental matrices
for i=1:n
    P=rand(3,4);
    while(rank(P)<3)
        P=rand(3,4);
    end
    Cams(:,:,i)=P/norm(P,'fro')*sc;
end

% generate noiseless fundamental matrices from cameras
Funds=zeros(3,3,m);
for k=1:m
    i=E(k,1); % left node
    j=E(k,2); % right node
    Pi=Cams(:,:,i); % left camera
    Pj=Cams(:,:,j); % right camera

    F=fundamental_det(Pi,Pj);
    Funds(:,:,k)=F/norm(F,'fro')*sc;
end

end