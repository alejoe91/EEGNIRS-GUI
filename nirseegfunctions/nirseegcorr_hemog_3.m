clear;
%% PARAMS
flag_print_fig = false;
numchan = 19;
figspath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
%% GET / PREPROCESS EEG DATA
% Input files
fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
fname = ['Omurtag327taskeyescloseboth_20140327_164238'];
badchans = {};
% PREPROCESS params
% data clip
dtbeg = 10; % sec
dtend = 10; % sec
% Bandpass filter and notch (normal bandpass range: 0.5-70 Hz)
filterlofreq = .5; % Hz
filterhifreq = 70; % Hz
filternotch = 60; % 0: no notch;
file1 = [fpath fname '.edf'];
[eegchanlabels, eegchanlocs, xxeegf, evteeg, Fs] = getpreprocess_eeg( file1, dtbeg, dtend, filterlofreq, filterhifreq, filternotch );
xxeegf = -xxeegf;
%% GET / PREPROCESS NIRS DATA
% Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
% Sample rate
Fsnirs = 6.25;
filterlofreq = .01; % Hz
filterhifreq = 1.1; % Hz
filternotch = 0; % 0: no notch;

for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end

fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
fpath2 = 'Hemoglobin Changes (NO FILTER)\Taskeyesclose\';

% TOO MANY CHANNELS HAVE BAD DATA
fname = 'NIRS-2014-03-27_004_deoxyhb_T1to2170_C1to20.txt';
% CZ PZ P4 BAD
fnameR = 'NIRS-2014-03-27_004_oxyhb_T1to2170_C1to20.txt';

badchans = {};

filenirs = [ fpath  fpath2  fname ];
filenirsR = [ fpath  fpath2  fnameR ];

% LOAD DATA
mydat = load( [ fpath  fpath2  fname  ]  );
mydatR = load( [ fpath  fpath2  fnameR  ]  );

% OMIT UNUSED CHANS
%
% NOTES:
% THERE IS A CORRECTION HERE FROM 20-20 TO 17-17
% ONLY 1-1, 2-2, 3-3, etc.
%
[nsamples, nchans0] = size( mydat );
xxnirs = zeros( nsamples, 19 ); xxnirsR = zeros( nsamples, 19 );
for ichan = 1:nchans0,
    if 20==ichan, % we used 20-20 for chan 6
        xxnirs( : , 17 ) = mydat( : , ichan );
        xxnirsR( : , 17 ) = mydatR( : , ichan );
    else, % we used 1-1 for chan1, 2-2 for chan2 etc.
        xxnirs( : , ichan ) = mydat( : , ichan );
        xxnirsR( : , ichan ) = mydatR( : , ichan );
    end
end

% FILTER
xxnirsf = mybutterandiirnotchfilters( xxnirs, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
xxnirsfR = mybutterandiirnotchfilters( xxnirsR, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );

% NIRS CHAN LABELS
for ii = 1:19,
    nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
    fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
end

%% GET NIRS EVENT FILE
fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
fpath2 = 'Detector Readings\Taskeyesclose\';
fname = 'NIRS-2014-03-27_004';
filenirsevt = [ fpath  fpath2  fname  '.evt' ];
%% SYNC

% EEG event sample numbers
Vevthresh = 1000;
nevteeg = remove_contiguous(  find( Vevthresh<evteeg )  );
tteeg = [1:length(evteeg)]./Fs;
tevteeg = nevteeg./Fs;
figure(110);clf; plot( tteeg, evteeg ); hold on; plot( tevteeg, ones(length(nevteeg)).*Vevthresh,'r.'); %axis([60 65 -500 1200])
nbegeeg = nevteeg(1);

% NIRS event sample numbers
myevtdat = load( [ filenirsevt ]  );
nevtnirs = myevtdat(:,1);
ttnirs = [1:length(xxnirsf)]./Fsnirs;
tevtnirs = nevtnirs./Fsnirs;
figure(210);clf; plot( ttnirs, xxnirsf(:,1) ); hold on; plot( tevtnirs, ones(length(nevtnirs)).*.04,'r.')
nbegnirs = nevtnirs(1);

% slippage
tslip = tevteeg(1)-tevtnirs(1);
tevteeg
tevtnirs + tslip

% clip eeg data at first event
xxeegfclip = xxeegf( nbegeeg:end, : );
tteegclip = tteeg( nbegeeg:end );
tteegclip = tteegclip - tteegclip(1);

% clip nirs data at first event
xxnirsfclip = xxnirsf( nbegnirs:end, :);
xxnirsfRclip = xxnirsfR( nbegnirs:end, :);
ttnirsclip = ttnirs( nbegnirs:end );
ttnirsclip = ttnirsclip - ttnirsclip(1);

if ttnirsclip(end)>tteegclip(end),
    % since nirs is longer, clip end of nirs
    nnn = min( find( tteegclip(end)<ttnirsclip ) ) - 1;
    ttnirsclip = ttnirsclip(1:nnn); xxnirsfclip = xxnirsfclip(1:nnn, :);
else
    nnn = min( find( ttnirsclip(end)<tteegclip ) ) - 1;
    tteegclip = tteegclip(1:nnn); xxeegfclip = xxeegfclip(1:nnn, :);
end

%figure(330); clf; plot(tteegclip, xxeegfclip(:,1), ttnirsclip, -50 + xxnirsfclip(:, 1) )
%figure(332); clf;  plot(tteegclip, xxeegfclip(:,1), ttnirsclip, -50 + xxnirsfclip(:, 1), [ttnirsclip(nnn) ttnirsclip(nnn)], [-50 50],'r*' )


%% DOWNSAMPLE (INTERPOLATE) EEG TO NIRS DATA POINTS
for ii = 1:19,
    xxeegfclipds(:,ii) = spline( tteegclip, xxeegfclip(:,ii), ttnirsclip );
end
%figure(430); clf; plot(ttnirsclip, xxeegfclipds(:,1), ttnirsclip, -50 + xxnirsfclip(:,1) )
%% LAGGED CORR
Tw = 40; % analysis window (sec)
Nwmax = 100000000; % max num windows to use (large number implies use max available)
nlags = ceil( Fsnirs.*(Tw-1) ); % max number of lags to compute (samples)

numsampW = floor( Fsnirs.*Tw );
nstds = [];

maxdim = min( Nwmax, floor( length(xxnirsfclip)./numsampW ) );

% first window
iw = 1;
ibeg = (iw-1).*numsampW + 1;
iend = ibeg + numsampW - 1;

XCFmean = []; XCFRmean = [];

fprintf( 'Calculating lagged correlation maxima for all windows all channel pairs\n' );

% loop through windows
while iw<=Nwmax & iend<=length(xxnirsfclip),
    
    dat1 = xxeegfclipds( ibeg:iend, : );
    dat2 = xxnirsfclip( ibeg:iend, : );
    dat2R = xxnirsfRclip( ibeg:iend, : );
    
    dat1 = mynormalize( dat1 );
    dat2 = mynormalize( dat2 );
    dat2R = mynormalize( dat2R );
    
    for ii = 1:numchan,
        chan1 = upper( nirschanlabels(ii) );
        for jj = 1:numchan,
            chan2 = upper ( eegchanlabels(jj) );
            
            if strcmp( chan1, chan2 ),
                fprintf( 'iw %5d/%5d  %6s %6s \n', iw, maxdim, char(chan1), char(chan2) );
                
                [ XCF, Lags ] = xcorr(  dat1(:,ii), dat2(:,jj),  nlags, 'coeff' );
                [ XCFR, Lags ] = xcorr(  dat1(:,ii), dat2R(:,jj),  nlags, 'coeff' );
                
                XCFmean = [ XCFmean, XCF ];
                XCFRmean = [ XCFRmean, XCFR ];
                
                if false, % plot correlation vs. lag
                    figure(880); clf;
                    plot( Lags./Fsnirs, XCF, 'r' ); hold on;
                    plot( Lags./Fsnirs, XCFR, 'g' ); hold on;
                    maxlim = max(Lags./Fsnirs);
                    axis([-maxlim  maxlim  -.5  .5]);
                    grid on;
                    titlestr = [  'iw ' num2str(iw) ', chans ' char(chan1) ' ' char(chan2) ];
                    title( titlestr );
                    pause(1);
                end
                
            end
        end
    end
    iw = iw + 1;
    % bounds of data for next window
    ibeg = (iw-1).*numsampW + 1;
    iend = ibeg + numsampW - 1;
end
Nw = iw - 1; % Total num of windows

%
maxlim = max(Lags./Fsnirs);
ymaxlim = .4;
XCFmean = mean( XCFmean, 2 );
XCFRmean = mean( XCFRmean, 2 );

%legend( [ pp1 pp2 pp1d pp2d],'eeg-eeg','nirs-nirs','eeg-laggedeeg','nirs-laggednirs'); legend boxoff

maxlim = max(Lags./Fsnirs);

%%
fontsz = 22; fontsztic = 16; linsz = 3;
figure(990); clf;

% smooth
ns = 9;
XCFmeans = mysmooth(XCFmean,ns);
XCFRmeans = mysmooth(XCFRmean,ns);

xminlim = -0;
xmaxlim = 20;
yminlim = -.05;
ymaxlim = .05;

pp = plot( Lags./Fsnirs, XCFmeans, 'color', 'r' , 'linewidth', linsz); hold on;
pp = plot( Lags./Fsnirs, XCFRmeans,'color', 'b' , 'linewidth', linsz); hold on;
axis(  [xminlim  xmaxlim  yminlim  ymaxlim]  );

xlabel('Time Lag (s)', 'fontsize', fontsz); ylabel('Correlation', 'fontsize', fontsz);
set( gca, 'fontsize', fontsztic );

set(gca,'xTick',[xminlim:10:xmaxlim]);
set(gca,'yTick',[yminlim:.025:ymaxlim]);

%%
if flag_print_fig,
    set( gcf, 'PaperPositionMode', 'auto')
    print(  '-dtiff', '-zbuffer', '-r600',  [figspath  'fig_task_eeg_nirs_laggedcorr'] );
end

%%
%{
TODO
- migrate the other progs fseeg fsnirs to the above subroutines
%}






