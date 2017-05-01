%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot gmm model pdf versus data histogram
%
function ok=plot_gmm(data,mu_gmm,sig_gmm,pp_gmm)

KS=length(ww_gmm);
fit_pdf=0*data;
for kks=1:KS
    fit_pdf=fit_pdf+pp_gmm(kks)*normpdf(data,mu_gmm(kks),sig_gmm(kks));
    plot(data,pp_gmm(kks)*normpdf(data,mu_gmm(kks),sig_gmm(kks)),'g');
end

plot(data,fit_pdf,'r');

grid on; 

ylabel('hist/pdf');
xlabel('x'); 

ok=1;

