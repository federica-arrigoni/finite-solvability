%% Compute the fundamental matrix given two cameras

function F=fundamental_det(P,Q)

F=zeros(3);
for h=1:3
    for k=1:3
        Pk=P; Pk(k,:)=[];
        Qh=Q; Qh(h,:)=[];
        F(h,k)=(-1)^(h+k)*det([Pk;Qh]);
    end
end

end

