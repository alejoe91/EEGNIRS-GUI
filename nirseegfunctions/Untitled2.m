


clear;
Fs = 250;
freq = .5;
tt = [1:2000]./Fs;
xx = sin(2.*pi.*freq.*tt) + .1.*randn(1,length(tt));
figure(110);clf;
plot(tt,xx); hold on;
figure(124);clf;
yy = mysmooth( xx,9);
plot(tt, yy,'b');


