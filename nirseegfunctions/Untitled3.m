clear;

Fs = 250;
freq = .5;
tt = [1:2000]./Fs;
xx = sin(2.*pi.*freq.*tt+pi./2) + .1.*randn(1,length(tt));
yy = sin(2.*pi.*freq.*tt) + .1.*randn(1,length(tt));

[cc,lags] = xcorr(xx,yy, 4.*Fs,'coeff');

figure(110);clf;
plot(tt,xx); hold on;
plot(tt,yy,'r');

figure(120);clf;
plot(lags,cc);
