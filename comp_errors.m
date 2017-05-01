%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for computing errors
%
function C_ERR=comp_errors(mi_true,pp_true,sig_true,mi_est,N_samp)

g_error=0;
KS=length(mi_true);
for kkk=1:KS
      g_error=g_error+min(abs(mi_true(kkk)-mi_est))*(sqrt(N_samp*pp_true(kkk))/(sig_true(kkk)));
end
C_ERR=g_error/KS;
