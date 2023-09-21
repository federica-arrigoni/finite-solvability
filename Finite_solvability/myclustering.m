%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023

function  [fcomp,C] = myclustering(U, tol)
% Cluster together rows of U that are equal up to tol
% Each row corresponds to one edge in the viewing graph
% Simple clustering with fixed radius "tol": candidate centers are all the
% rows, starting from the 1st and proceeding in a greedy fashion

remaining = 1:size(U,1);
fcomp = nan(size(U,1),1);
ind_comp = 1;

while (~isempty(remaining))

    % centre of the component: first out of edges not assigned
    C(ind_comp,:)= U(remaining(1),:); 
    
    % find all the rows that are equal to the current one
    edges_in_component = vecnorm(U(remaining(1),:) - U, 2, 2) < tol ;
    fcomp( edges_in_component ) = ind_comp;
    
    remaining = find(isnan(fcomp)); % update edges not assigned to any component
    ind_comp = ind_comp +1; % go to the next component
    
end

assert(all(isfinite(fcomp)));

end