function plot7maps_horiz( fignum, titlestrIN, Evec, Lam, chanstoinclude, labels, xchan );

fh = figure(fignum);clf;
set(fh,'color','w');

iverbose = 0; labelclr = 'k';
ivecnum = 1;
for jjj = 1:7,
  %  for kkk = 1:2,
        s2 = subplot(1,7,ivecnum);
        rendermontage3( 1234, labels, xchan, iverbose, labelclr, chanstoinclude, Evec(:,ivecnum), .0 );
        titlestr0 = [ num2str(ivecnum) ': ' num2str(Lam(ivecnum)) ];
        if 1==ivecnum, titlestr0 = [titlestrIN ' mode ' titlestr0]; end
        
        title(titlestr0,'fontsize',12);
        
        s2Pos = get(s2,'position');
        ivecnum = ivecnum + 1;
 %   end
end
% colorbar
hb = colorbar('location','southoutside');
set(s2,'position',s2Pos);
%set(hb,'XTickLabel',-1:.15:1);
%set(hb,'XTick',[]);
%xlim(hb, [-1 1]);
caxis([-1, 1]);
%get(hb);
set(hb,'Box','off');
set(hb,'LineWidth',.01);
%set(hb,'XColor','w');
V = get(hb,'OuterPosition');
set(hb,'OuterPosition',[V(1) V(2)./2 V(3) V(4)./1.5]);
U = get(hb,'TickLength');
set(hb,'TickLength',[U(1)./9 U(2)./9]);
set(hb,'TickDir','out');