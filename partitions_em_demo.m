% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%       partitions em demo
%
% script file for demonstratig comparison of algorithms for comparisons 
% of initialization methods for EM algorithm
            
% Mixture paramters
%
% number of components
KS=10;
%
% draw standard deviations of true mixture components
%sig_v=(1:KS)/KS;
sig_true=unifrnd(0.05,1,1,KS);
sig_true=sig_true(randperm(KS));
%
% component wieghts
pp_true=1:KS;
pp_true=pp_true/sum(pp_true);
pp_true=pp_true(randperm(KS));
%
% overlap coefficient
ov=0.15;
%
% compute mixture components expectations by using assumed overlaps
mu_true=zeros(1,KS);
for kk=2:KS
        mu_true(kk)=mu_true(kk-1)+(-2*log(ov))*sqrt(sig_true(kk-1)^2+sig_true(kk)^2);
end

% generate mixture sample with N elemnts
N=1000;
data=g_mix_gen(mu_true,sig_true,pp_true,N);
data=sort(data)';
      
            
% data buffers for initial values of parameters
%
% method - equal quantilles (inverse CDF) - EQ
sig_ini_EQ=zeros(1,KS);
pp_ini_EQ=zeros(1,KS);
mu_ini_EQ=zeros(1,KS);
part_sticks_EQ=zeros(1,KS-1);
%
% method - hierarchical clustering average linkage - h_clu_a
sig_ini_hclu_a=zeros(1,KS);
pp_ini_hclu_a=zeros(1,KS);
mu_ini_hclu_a=zeros(1,KS);
part_sticks_hclu_a=zeros(1,KS-1);
%
% method - hierarchical clustering complete linkage - h_clu_c
sig_ini_hclu_c=zeros(1,KS);
pp_ini_hclu_c=zeros(1,KS);
mu_ini_hclu_c=zeros(1,KS);
part_sticks_hclu_c=zeros(1,KS-1);   
%
% method - dynamic programmin, version Q4
sig_ini_dp_4=zeros(1,KS);
pp_ini_dp_4=zeros(1,KS);
mu_ini_dp_4=zeros(1,KS);
part_sticks_dp_4=zeros(1,KS-1);   

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute initial condition by using equal quantilles (invCDF)
% draw partition
% estimate parameters - draw estimated mixture pdf versus histogram
% compute log likelihood and errors
for kkp=1:KS
    pp_ini_EQ(kkp)=1/KS;
    mu_ini_EQ(kkp)=mean(data(round((kkp-1)*N/KS)+1:round((kkp)*N/KS)));   
    sig_ini_EQ(kkp)=std(data(round((kkp-1)*N/KS)+1:round((kkp)*N/KS)));  
end
% store partition
for kkp=1:KS-1
    part_sticks_EQ(kkp)=data(round((kkp)*N/KS));
end
%
% draw obtained partition
figure(1);
subplot(4,1,1);
ok=draw_part(data,part_sticks_EQ);
title('Partitions by equal quantilles - EQ')
%
% estimate mixture parameters and store results
[mu_est_EQ,sig_EQ_cdf,pp_est_EQ,l_lik_EQ] = g_mix_est_fast_lik(data,KS,mu_ini_EQ,sig_ini_EQ,pp_ini_EQ);
%
% draw data histogram versus fitted mixture model
figure(2);
subplot(4,1,1);
ok=draw_hist_pdf(data,mu_est_EQ,sig_EQ_cdf,pp_est_EQ);
% compute error between true and estimated parameters
C_ERR_EQ=comp_errors(mu_true,pp_true,sig_true,mu_est_EQ,N);
title(['Histogram versus estimated pdf - EQ,  l-lik= ' num2str(l_lik_EQ)  '  D=  ' num2str(C_ERR_EQ)])
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute initial condition by using hclu_a
% draw partition
% estimate parameters - draw estimated mixture pdf versus histogram
clusters=h_clu_a(data,KS);
for kkp=1:KS
    pp_ini_hclu_a(kkp)=(clusters(kkp,2)-clusters(kkp,1))/N;
    mu_ini_hclu_a(kkp)=mean(data(clusters(kkp,1):clusters(kkp,2)));   
    sig_ini_hclu_a(kkp)=std(data(clusters(kkp,1):clusters(kkp,2)));  
end
%
% store partition
for kkp=1:KS-1
    part_sticks_hclu_a(kkp)=data(clusters(kkp,2));
end
%
% draw obtained partition
figure(1);
subplot(4,1,2);
ok=draw_part(data,part_sticks_hclu_a);
title('Partitions by hierarchical clust. average - h-clu-a')
%
% estimate mixture parameters and store results
[mu_est_hclu_a,sig_est_hclu_a,pp_est_hclu_a,l_lik_hclu_a] = g_mix_est_fast_lik(data,KS,mu_ini_hclu_a,sig_ini_hclu_a,pp_ini_hclu_a);
%
% draw data histogram versus fitted mixture model
figure(2);
subplot(4,1,2);
ok=draw_hist_pdf(data,mu_est_hclu_a,sig_est_hclu_a,pp_est_hclu_a);
% compute error between true and estimated parameters
C_ERR_hclu_a=comp_errors(mu_true,pp_true,sig_true,mu_est_hclu_a,N);
title(['Histogram versus estimated pdf - h-clu-a,  l-lik= ' num2str(l_lik_hclu_a) '  D=  ' num2str(C_ERR_hclu_a)])


            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute initial condition by using hclu_c
% draw partition
% estimate parameters - draw estimated mixture pdf versus histogram
clusters=h_clu_c(data,KS);
for kkp=1:KS
    pp_ini_hclu_c(kkp)=(clusters(kkp,2)-clusters(kkp,1))/N;
    mu_ini_hclu_c(kkp)=mean(data(clusters(kkp,1):clusters(kkp,2)));   
    sig_ini_hclu_c(kkp)=std(data(clusters(kkp,1):clusters(kkp,2)));  
end
%
% store partition
for kkp=1:KS-1
    part_sticks_hclu_c(kkp)=data(clusters(kkp,2));
end
%
% draw obtained partition
figure(1);
subplot(4,1,3);
ok=draw_part(data,part_sticks_hclu_c);
title('Partitions by hierarchical clust. complete - h-clu-c')
%
% estimate mixture parameters and store results
[mu_est_hclu_c,sig_est_hclu_c,pp_est_hclu_c,l_lik_hclu_c] = g_mix_est_fast_lik(data,KS,mu_ini_hclu_c,sig_ini_hclu_c,pp_ini_hclu_c);
%
% draw data histogram versus fitted mixture model
figure(2);
subplot(4,1,3);
ok=draw_hist_pdf(data,mu_est_hclu_c,sig_est_hclu_c,pp_est_hclu_c);
% compute error between true and estimated parameters
C_ERR_hclu_c=comp_errors(mu_true,pp_true,sig_true,mu_est_hclu_c,N);
title(['Histogram versus estimated pdf - h-clu-c,  l-lik= ' num2str(l_lik_hclu_a) '  D=  ' num2str(C_ERR_hclu_c)])

             
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute initial condition by using dynpro - fourth version, Q4
% 
% 
ver=4;
[Q,opt_part]=dyn_pr_split(data,KS-1,ver);
part_cl=[1 opt_part N+1];             
for kkp=1:KS
    pp_ini_dp_4(kkp)=(part_cl(kkp+1)-part_cl(kkp))/N;
    mu_ini_dp_4(kkp)=mean(data(part_cl(kkp):part_cl(kkp+1)-1));   
    sig_ini_dp_4(kkp)=std(data(part_cl(kkp):part_cl(kkp+1)-1));  
end
%
% store partition
for kkp=1:KS-1
    part_sticks_dp_4(kkp)=data(part_cl(kkp+1));
end
%
% draw obtained partition
figure(1);
subplot(4,1,4);
ok=draw_part(data,part_sticks_dp_4);
title('Partitions by dynamic programming - dp-4')
%
% estimate mixture parameters and store results
[mu_est_dp_4,sig_est_dp_4,pp_est_dp_4,l_lik_dp_4] = g_mix_est_fast_lik(data,KS,mu_ini_dp_4,sig_ini_dp_4,pp_ini_dp_4);
%
% draw data histogram versus fitted mixture model
figure(2);
subplot(4,1,4);
ok=draw_hist_pdf(data,mu_est_dp_4,sig_est_dp_4,pp_est_dp_4);
% compute error between true and estimated parameters
C_ERR_dp_4=comp_errors(mu_true,pp_true,sig_true,mu_est_dp_4,N);
title(['Histogram versus estimated pdf - dp-4,  l-lik= ' num2str(l_lik_dp_4)  '  D=  ' num2str(C_ERR_dp_4)])
          
