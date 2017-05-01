%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute sample partition by dynamic programming 
%
function [Q,opt_part]=dyn_pr_split(data,K_gr,ver)

% initialize
Q=zeros(1,K_gr);
N=length(data);
p_opt_idx=zeros(1,N);
p_aux=zeros(1,N);
opt_pals=zeros(K_gr,N);
for kk=1:N;
    p_opt_idx(kk)=my_qu_ix(data(kk:N),ver);
end


% aux_mx
aux_mx=zeros(N,N);
for kk=1:N-1
   for jj=kk+1:N
       aux_mx(kk,jj)= my_qu_ix(data(kk:jj-1),ver);
   end
end

% iterate
for kster=1:K_gr
   % kster
   for kk=1:N-kster
       for jj=kk+1:N-kster+1
           p_aux(jj)= aux_mx(kk,jj)+p_opt_idx(jj);
       end
       [mm,ix]=min(p_aux(kk+1:N-kster+1));
       p_opt_idx(kk)=mm;
       opt_pals(kster,kk)=kk+ix(1);
   end
   Q(kster)=p_opt_idx(1);
end


% restore optimal decisions
opt_part=zeros(1,K_gr);
opt_part(1)=opt_pals(K_gr,1);
for kster=K_gr-1:-1:1
   opt_part(K_gr-kster+1)=opt_pals(kster,opt_part(K_gr-kster));
end



% auxiliary function
function wyn=my_qu_ix(invec,ver)
switch ver
    case 1
         if length(invec)<3
           wyn=inf;
        else
           wyn=var(invec);
        end 
    case 2
        if length(invec)<3
           wyn=inf;
        else
            wyn=std(invec);
        end
    case 3
        par=0.0;
        if length(invec)<3
           wyn=inf;
        else
           wyn=(par+std(invec))/(max(invec)-min(invec));
        end
    case 4
        par=0.1;
        if length(invec)<3
           wyn=inf;
        else
           wyn=(par+std(invec))/(max(invec)-min(invec));
        end
end


