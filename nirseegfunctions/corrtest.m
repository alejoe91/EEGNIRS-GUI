clear;format compact;
%% orig signal
T = 2000;
Fs = 6.25;
N = floor( T.*Fs );
dt = 1./Fs;
tt = [1:N]./Fs;
anoise_orig = 2.8;
anoise_delay = .3;
anoise_conv = .1;
if 0,
    % ornstein uhlenbeck process
    tau = .25; xx = zeros(1,N); for ii = 2:N, xx(ii) = xx(ii-1) + dt.*(-1./tau.*xx(ii-1) + randn);end
elseif 1,
    % sinusoid
    freq = 1./75; % Hz
    xx = sin( 2.*pi.*freq.*tt );
else,
    % random
    xx = randn(1,N);
    %xx(20:30) = xx(20:30) + 8;
end
%% lagged corr params
Tw = 120;
nw = ceil( Tw.*Fs );
twindowdisp = 60;
%% add noise to orig signal
xx = xx + anoise_orig .* randn( size(xx) );
%% delayed sig
tdelay = 4;
ndelay = floor(tdelay.*Fs);
xx_del = xx;
xx_del(ndelay:end) = xx(1:end-ndelay+1);
xx_del = xx_del + anoise_delay .* randn( size(xx_del) );
figure(100); clf; plot(tt,xx,'r',tt,xx_del,'g');  legend(['orig'],'delayed'); legend boxoff
axis([0 twindowdisp min([xx]) max(xx)])
%% impulse resp
tconv = 20;
M = 2.*ceil(tconv.*Fs) + 1;
ttconv = [1:M]./Fs - tconv - dt;
tmean = 9; %nmean = floor(tmean.*Fs);
tstd = tmean./2; %nstd
fff = zeros(M,1);
fff = exp( -(ttconv - tmean ).^2 ./(2.*tstd) );
%fff = fff./sum(fff);
figure(200);clf; plot( ttconv,fff);title('Impulse response');
%% convovled signal
xx_conv = conv(xx, fff, 'same');
xx_conv = xx_conv + anoise_conv .* randn( size(xx_conv) );
figure(300); clf; plot(tt,xx,'r',tt,xx_conv,'g','linewidth',2);  legend('orig','convolved'); legend boxoff; axis([0 twindowdisp min([xx xx_conv]) max([xx xx_conv])]);
figure(320); clf; plot(tt,xx_conv,'g','linewidth',2);  legend('convolved'); legend boxoff; axis([20 20+twindowdisp min([ xx_conv]) max([ xx_conv])]);
%% lagged corr
[corrs_del, lags] = xcorr( xx_del, xx, nw, 'coeff' );
[corrs_conv, lags] = xcorr( xx_conv, xx, nw, 'coeff' );
[corrs_auto_orig, lags] = xcorr( xx, xx, nw, 'coeff' );
figure(400);clf;
pp1=plot( lags./Fs, corrs_del,'b'); hold on; title(['Lagged Correlation, delay ' num2str(tdelay)] );
pp2=plot( lags./Fs, corrs_conv,'r');grid on;
pp3=plot( lags./Fs, corrs_auto_orig,'k');grid on;
legend([pp1 pp2 pp3],'Delayed','Convolved','Auto (orig)'); legend boxoff
%% lagged mutual info
numlagsamples = nw;
[mi_del, milags_del] = laggedmutualinfo( xx, xx_del, numlagsamples );
[mi_conv, milags_conv] = laggedmutualinfo( xx, xx_conv, numlagsamples );
milags_del = milags_del./Fs;
milags_conv = milags_conv./Fs;
figure(500);clf;
subplot(2,1,1), plot( milags_del, mi_del, 'b');hold on;title(['Lagged Mutual Information - Delayed ' num2str(tdelay)] );hold on;grid on;
subplot(2,1,2), plot( milags_conv, mi_conv, 'r');title(['Lagged Mutual Information - Convolved ' num2str(tmean)] );
hold on;grid on;






