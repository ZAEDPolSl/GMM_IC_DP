%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% draw partition
%
function ok=draw_part(data,part_sticks)
[hn,xx]=hist(data,100);
N=length(data);
KSS=length(part_sticks);
shn=hn/(N*(xx(2)-xx(1)));
mxh=max(shn);
hold off
bar(xx,shn,1,'w');
hold on
for kkp=1:KSS
   plot([part_sticks(kkp) part_sticks(kkp)], [0 mxh], 'r');   
end
ok=1;