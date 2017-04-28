clear;
Fs = 250;
TT = 20;
f1 = .1;
f2 = 10;


N = TT.*Fs;
tt = [1:N]./Fs;

aa = .3;
A = sin(2.*pi.*f1.*tt);
xx = A.* sin(2.*pi.*f2.*tt) + aa.*randn(size(tt));

figure(220);clf; plot(tt,xx);


nw = 100;
num = nw-1;
denom = nw;
noverlap = floor( nw.* (num./denom) );
F = [];
Fs = 250;
[S,F,T,P] = spectrogram(xx,nw,noverlap,F,Fs);

indxs = find( 0<=F & F<=4 ); p_d = sum( P(indxs,:), 1 );
indxs = find( 8<=F & F<=12 ); p_a = sum( P(indxs,:), 1 );
indxs = find( 5<=F & F<=8 ); p_t = sum( P(indxs,:), 1 );
indxs = find( 13<=F & F<=18 ); p_b = sum( P(indxs,:), 1 );

figure(330);clf; plot(T,p_a,T,p_d,T,p_t,T,p_b);

figure(110);clf;
zz = 10.*log10(P);
maxz = max(max(zz));
minz = min(min(zz));
surf(T,F,zz,'edgecolor','none'); 
axis([0 TT 0 60 minz maxz] );
%axis tight
view(0,90)
xlabel('Time (s)'); ylabel('Hz');



