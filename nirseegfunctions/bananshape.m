clear;

%%
L = -20;
R = 40;
xmin = L;
xmax = R;
ymin = L;
ymax = R;
zmin = L;
zmax = R;

N = 100;
xxx = linspace(xmin,xmax,N);
yyy = linspace(ymin,ymax,N);
zzz = linspace(zmin,zmax,N);

k = .23;
d = 30;

for ii = 1:N,
    for jj = 1:N,
        for kk = 1:N,
            xx = xxx(ii); yy = yyy(jj); zz = zzz(kk);
            if zz<0,
                numer = zz.^2 .* exp( -k.*( (xx.^2+yy.^2+zz.^2).^.5  + ((xx-d).^2+yy.^2+zz.^2).^.5 ) );
                denom = (xx.^2+yy.^2+zz.^2).^(3./2) .* ...
                    ( ((xx-d).^2+yy.^2+zz.^2).^(3./2) ) .* ...
                    (k.* ((xx.^2+yy.^2+zz.^2).^.5+1)).*(k.* ( (xx-d).^2+yy.^2+zz.^2).^.5 + 1);
                P(ii,jj,kk) = numer./denom;
            else,
                P(ii,jj,kk) = 0;
            end
        end
    end
end
%%
[X,Y,Z] = meshgrid(xxx,yyy,zzz);

pcut = 8e-11;

figure(210);clf;
p = patch(isosurface(X,Y,Z,P,pcut));
isonormals(X,Y,Z,P,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
%view(3);
%axis tight
camlight
%lighting gouraud
grid on
xlabel('x'); ylabel('y'); zlabel('z')
axis([ xmin xmax ymin ymax zmin 0] )
view(90,0)
%%
Pflat = reshape(P,N^3,1,1);
[phist,xhist]=hist(Pflat,80);
figure(10);clf;
loglog(xhist,phist,'.-');

return;
