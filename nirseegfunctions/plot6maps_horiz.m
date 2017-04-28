function plot6maps_horiz( fignum, titlestrIN, Evec, Lam, chanstoinclude, labels, xchan, labelclr );

fh = figure(fignum);clf;
set(fh,'color','w');

iverbose = 0; 
ivecnum = 1;
for jjj = 1:6,
  %  for kkk = 1:2,
        s2 = subplot(1,6,ivecnum);
        if 6==jjj, displabel = 1; else, displabel = 0; end
        cutoff = 0;
        rendermontage3( 1234, labels, xchan, iverbose, labelclr, chanstoinclude, Evec(:,ivecnum), cutoff, displabel );
        titlestr0 = [ num2str(ivecnum) ': ' num2str(Lam(ivecnum)) ];
        if 1==ivecnum, titlestr0 = [titlestrIN ' mode ' titlestr0]; end
        
        title(titlestr0,'fontsize',12);
        
        s2Pos = get(s2,'position');
        ivecnum = ivecnum + 1;
 %   end
end
% colorbar
hb = colorbar('location','eastoutside');
set(s2,'position',s2Pos);
%set(hb,'XTickLabel',-1:.15:1);
%set(hb,'XTick',[]);
%xlim(hb, [-1 1]);
caxis([-1, 1]);
%get(hb);
set(hb,'Box','off');
set(hb,'LineWidth',.01);
%set(hb,'XColor','w');
U = get(hb,'TickLength');
set(hb,'TickLength',[U(1)./9 U(2)./9]);
set(hb,'TickDir','out');
V = get(hb,'OuterPosition');
set(hb,'OuterPosition',[V(1).*1.06 V(2) V(3)./1.1 V(4)./1]);


