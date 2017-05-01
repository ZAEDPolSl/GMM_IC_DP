%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sample partition by complete linkage hierarchical clustering
% 
function clusters=h_clu_c(sample,K_cl)

N=length(sample);

% actual number of clusters
N_act=N;

% distance vector
dist_vec=sample(2:N)-sample(1:N-1);

clusters=[(1:N)' (1:N)'];


% loop steps
for kk=1:N-K_cl

   % one step
   [minnn,ixmin]=min(dist_vec);
   N_act=N_act-1;

   ix_merge=ixmin(1);

   clusters(ix_merge,2)=clusters(ix_merge+1,2);
   clusters(ix_merge+1:N_act,:)=clusters(ix_merge+2:N_act+1,:);
   dist_vec(ix_merge:N_act-1)=dist_vec(ix_merge+1:N_act);

   
    if ix_merge > 1
       dist_vec(ix_merge-1)=udist_c(sample(clusters(ix_merge-1,1):clusters(ix_merge-1,2)),sample(clusters(ix_merge,1):clusters(ix_merge,2)));
    end
    if ix_merge < N_act
       dist_vec(ix_merge)=udist_c(sample(clusters(ix_merge,1):clusters(ix_merge,2)),sample(clusters(ix_merge+1,1):clusters(ix_merge+1,2)));
    end

    dist_vec=dist_vec(1:N_act-1); 
    clusters=clusters(1:N_act,:);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function dist=udist_c(kla1,kla2)
% complete linkage
dist=max(kla2)-min(kla1);

