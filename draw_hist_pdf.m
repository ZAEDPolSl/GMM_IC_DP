function ok=draw_hist_pdf(data,mu_gmm,sig_gmm,pp_gmm)
% draw histogram versus pdf

[hn,xx]=hist(data,100);
N=length(data);
shn=hn/(N*(xx(2)-xx(1)));
hold off
bar(xx,shn,1,'w');
hold on

KS=length(pp_gmm);
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

