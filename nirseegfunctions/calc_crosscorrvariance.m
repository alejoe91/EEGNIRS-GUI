function [ss] = calc_crosscorrvariance( x, y, maxlags, varlags )
%{
kramer, eden, cash, et al. 2009, physical review
eq. 2

x, y: two signals presumably correlated with some noise added
maxlags: the total lag to use in computing lagged correlations (samples)
varlags: the total lag to use in computing variance estimate
%}
[ cxx, Lags ] = xcorr(  x,  maxlags, 'coeff' );
[ cyy, Lags ] = xcorr(  y,  maxlags, 'coeff' );
[ cyx, Lags ] = xcorr(  x, y,  maxlags, 'coeff' );

if varlags>=maxlags, error(  ' **** varlags cannot be greater than maxlags' ); end;

ss = sum(cxx.*cyy)./(maxlags - varlags);

fsnirs = 6.25;
figure(9010); clf; 
plot(Lags./fsnirs,cxx, Lags./fsnirs,cyy ,Lags./fsnirs,cyx); 
hold on; plot(varlags./fsnirs,ss,'r');


fsdafads=987789;





