clear;
%% PARAMS
numchan = 19;
figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
%% GET prepocessed EEG DATA
% Input files
fpath = 'C:\downloads\nirseegdata\';
% fname = 'AhmetOmurtag-eyesclose1_20131216_193358.edf';
%  fname = 'AhmetOmurtag-eyesopen1_20131216_191542.edf';
%  fname = 'AhmetOmurtag-closenotrigger_20131216_202346.edf';
% fname = 'AhmetOmurtag-eyesopen3_20131216_194123.edf';
  fname = 'AhmetOmurtag-taskvft_20131216_200735.edf';
% PREPROCESS params
% data clip
dtbeg = 5; % sec
dtend = 5; % sec
% Bandpass filter and notch (normal bandpass range: 0.5-70 Hz)
filterlofreq = .5; % Hz
filterhifreq = 70; % Hz
filternotch = 60; % 0: no notch;
file1 = [fpath fname];
[eegchanlabels, eegchanlocs, xxeegf, evteeg, Fs] = getpreprocess_eeg( file1, dtbeg, dtend, filterlofreq, filterhifreq, filternotch );

%% GET PREPROCESSED NIRS DATA
% Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
% Sample rate
Fsnirs = 3.125;
filterlofreq = .01; % Hz
filterhifreq = .8; % Hz
filternotch = 0; % 0: no notch;

fpath = 'C:\downloads\nirseegdata\LabData\LabData\AhmetOmurtag1216\NIRS\';

if 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_001 - eyesopen1\';
    fname = 'NIRS-2013-12-16_001'; 
elseif 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_002 - eyesopen2\';
    fname = 'NIRS-2013-12-16_002';
elseif 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_004 - eyesopen3\';
    fname = 'NIRS-2013-12-16_004';
elseif 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_003 - eyesclose1\';
    fname = 'NIRS-2013-12-16_003'; 
elseif 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_005 - eyesclose2\';
    fname = 'NIRS-2013-12-16_005';
elseif 1
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_006 - taskvft\';
    fname = 'NIRS-2013-12-16_006'; 
elseif 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_007 - opennotrigger\';
    fname = 'NIRS-2013-12-16_007';
elseif 0
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_008 - closenotrigger\';
    fname = 'NIRS-2013-12-16_008';
end

fnamext = 'wl2';

filenirs = [ fpath  fpath2  fname  '.'  fnamext];
filenirsevt = [ fpath  fpath2  fname  '.evt' ];

[nirschanlabels, xxnirsf ] = getpreprocess_nirs(filenirs,filterlofreq,filterhifreq, filternotch, Fsnirs);

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

figure(330); clf; plot(tteegclip, xxeegfclip(:,1), ttnirsclip, -50 + xxnirsfclip(:, 1) )
%figure(332); clf;  plot(tteegclip, xxeegfclip(:,1), ttnirsclip, -50 + xxnirsfclip(:, 1), [ttnirsclip(nnn) ttnirsclip(nnn)], [-50 50],'r*' )


%% DOWNSAMPLE (INTERPOLATE) EEG TO NIRS DATA POINTS
for ii = 1:19,
    xxeegfclipds(:,ii) = spline( tteegclip, xxeegfclip(:,ii), ttnirsclip );
end
figure(430); clf; plot(ttnirsclip, xxeegfclipds(:,1), ttnirsclip, -50 + xxnirsfclip(:,1) )
%% LAGGED CORR
Tw = 250; % analysis window (sec)
Nwmax = 100000000; % max num windows to use (large number implies use max available)
nlags = ceil( Fsnirs.*30 ); % max number of lags to compute (samples)

numsampW = floor( Fsnirs.*Tw );
nstds = [];

maxdim = min( Nwmax, floor( length(xxnirsfclip)./numsampW ) );

% max correlation
corrmax = zeros( numchan, numchan, maxdim );
% lag at which corr is maximum
lagcorrmax = zeros( numchan, numchan, maxdim );

% first window
iw = 1;
ibeg = (iw-1).*numsampW + 1;
iend = ibeg + numsampW - 1;

XCFmean = []; XCFmeansurr = [];

fprintf( 'Calculating lagged correlation maxima for all windows all channel pairs\n' );

% loop through windows
while iw<=Nwmax & iend<=length(xxnirsfclip),

    dat1 = xxeegfclipds( ibeg:iend, : ); 
    dat2 = xxnirsfclip( ibeg:iend, : );
    
    dat1 = mynormalize( dat1 ); dat2 = mynormalize( dat2 );
    
    Tadv = 3; % sec
    nlagsurrogate = ceil( Fsnirs.*Tadv ); % lag of surrogate data (samples)
    fprintf('Surrogate data lag: %f seconds\n', nlagsurrogate./Fsnirs);
%    dat2surrogate = [ dat1( nlagsurrogate:end, : ); dat1( 1:nlagsurrogate-1, : ) ];      % advanced by Tadv secs
    dat1surrogate = [ dat1( length(dat1)-nlagsurrogate+1:end, : ); dat1( 1:length(dat1)-nlagsurrogate , : ) ]; % lagged by Tadv secs
    dat2surrogate = [ dat2( length(dat2)-nlagsurrogate+1:end, : ); dat2( 1:length(dat2)-nlagsurrogate , : ) ]; % lagged by Tadv secs
        
    for ii = 1:numchan,
        chan1 = upper( nirschanlabels(ii) );
        for jj = 1:numchan,
            chan2 = upper ( eegchanlabels(jj) );

            if strcmp( chan1, chan2 ),
                fprintf( 'iw %5d/%5d  %6s %6s \n', iw, maxdim, char(chan1), char(chan2) );

                [ XCF, Lags, Bounds ] = crosscorr(  dat1(:,ii), dat2(:,jj),  nlags, nstds );
                [cmx, indxmaxcorr] = max(XCF);
                corrmax(ii, jj, iw) = cmx;
                lagcorrmax(ii, jj, iw) = Lags(indxmaxcorr);
                
                % surrogate                
                [ XCFsurr1, Lags, Bounds ] = crosscorr(  dat1(:,ii), dat1(:,jj),  nlags, nstds );
                [ XCFsurr2, Lags, Bounds ] = crosscorr(  dat2(:,ii), dat2(:,jj),  nlags, nstds );
                [ XCFsurr1d, Lags, Bounds ] = crosscorr(  dat1(:,ii), dat1surrogate(:,jj),  nlags, nstds );
                [ XCFsurr2d, Lags, Bounds ] = crosscorr(  dat2(:,ii), dat2surrogate(:,jj),  nlags, nstds );
                
                XCFmean = [ XCFmean, XCF ];
                XCFmeansurr1 = [ XCFmeansurr, XCFsurr1 ];
                XCFmeansurr2 = [ XCFmeansurr, XCFsurr2 ];
                XCFmeansurr1d = [ XCFmeansurr, XCFsurr1d ];
                XCFmeansurr2d = [ XCFmeansurr, XCFsurr2d ];
                
                if false, % plot correlation vs. lag
                    figure(880); clf;
                    plot( Lags./Fsnirs, XCF, 'r' ); hold on;
                    plot( Lags./Fsnirs, XCFsurr, 'g' ); hold on;
                    maxlim = max(Lags./Fsnirs); 
                    axis([-maxlim  maxlim  -.5  .5]); 
                    grid on;
                    titlestr = [  'iw ' num2str(iw) ', chans ' char(chan1) ' ' char(chan2)  ',  maxcorr ' num2str( corrmax(ii, jj, iw) ) ', lag maxcorr: ' num2str(lagcorrmax(ii, jj, iw)  ) ];
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
XCFmean = mean( XCFmean, 2 ); XCFmeansurr = mean( XCFmeansurr, 2 );
figure(890); clf;
subplot(2,2,1),
pp1 = plot( Lags./Fsnirs, XCFmeansurr1, 'g' ); hold on;
pp = plot( Lags./Fsnirs, XCFmean, 'r' ); hold on; axis([-maxlim  maxlim  -ymaxlim  ymaxlim]); grid on;
title(['eeg-nirs and eeg-eeg']);
subplot(2,2,2),
pp2 = plot( Lags./Fsnirs, XCFmeansurr2, 'g' ); hold on;
pp = plot( Lags./Fsnirs, XCFmean, 'r' ); hold on; axis([-maxlim  maxlim  -ymaxlim  ymaxlim]); grid on;
title(['eeg-nirs and nirs-nirs']);

subplot(2,2,3)
pp1d = plot( Lags./Fsnirs, XCFmeansurr1d, 'g' ); hold on;
pp = plot( Lags./Fsnirs, XCFmean, 'r' ); hold on; axis([-maxlim  maxlim  -ymaxlim  ymaxlim]); grid on;
title(['eeg-nirs and eeg-eeglagged']);

subplot(2,2,4)
pp2d = plot( Lags./Fsnirs, XCFmeansurr2d, 'g' ); hold on;
pp = plot( Lags./Fsnirs, XCFmean, 'r' ); hold on; axis([-maxlim  maxlim  -ymaxlim  ymaxlim]); grid on;
title(['eeg-nirs and nirs-nirslagged']);

%legend( [ pp1 pp2 pp1d pp2d],'eeg-eeg','nirs-nirs','eeg-laggedeeg','nirs-laggednirs'); legend boxoff

maxlim = max(Lags./Fsnirs);

xlabel('Lag (s)'); ylabel('Correlation');

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






