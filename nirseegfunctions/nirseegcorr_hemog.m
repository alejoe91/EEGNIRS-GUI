clear;
%% PARAMS
flag_print_fig = false;
numchan = 19;
figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
%% GET prepocessed EEG DATA
% Input files
fpath = 'F:\LabData\LabData\AhmetOmurtag1216\EEG\';
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
xxeegf = -xxeegf;
%% GET PREPROCESSED NIRS DATA
% Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
% Sample rate
Fsnirs = 3.125;
filterlofreq = .01; % Hz
filterhifreq = .8; % Hz
filternotch = 0; % 0: no notch;

fpath = 'F:\LabData\LabData\AhmetOmurtag1216\NIRS\';

icycle = 6;

if 1==icycle,
    fpath2 = 'Hemoglobin Con.Changes (NOFilter) (NIRSLAB)\1-eyesopen\';
%    fpath2 = 'Hemoglobin Con. Changes  (NIRSLAB)\Ahmet Omurtag- 2013-12-16_001 - eyesopen1\';    
    fname = 'NIRS-2013-12-16_001_oxyhb_T1to1035_C1to20.txt';        
elseif 2==icycle,
    fpath2 = 'Hemoglobin Con.Changes (NOFilter) (NIRSLAB)\2-eyesopen2\';
    fname = 'NIRS-2013-12-16_002_deoxyhb_T1to1077_C1to20.txt';       
elseif 3==icycle,
    fpath2 = 'Hemoglobin Con.Changes (NOFilter) (NIRSLAB)\4-eyesopen3\';
    fname = 'NIRS-2013-12-16_004_deoxyhb_T1to1008_C1to20.txt';       
elseif 4==icycle,
    fpath2 = 'Hemoglobin Con.Changes (NOFilter) (NIRSLAB)\3-eyesclose1\';
    fname = 'NIRS-2013-12-16_003_oxyhb_T1to1077_C1to20.txt';         
    fnamefigconn = 'fig_nirs_conn_oxy_eyesclosed1';
elseif 5==icycle,
    fpath2 = 'Hemoglobin Con.Changes (NOFilter) (NIRSLAB)\5-eyesclose2\';
    fname = 'NIRS-2013-12-16_005_oxyhb_T1to1029_C1to20.txt';     
elseif 6==icycle,
    fpath2 = 'Hemoglobin Con.Changes (NOFilter) (NIRSLAB)\6-taskvft\';
    fname = 'NIRS-2013-12-16_006_oxyhb_T1to1042_C1to20.txt';
    fnamefigconn = 'fig_nirs_conn_oxy_taskvft';
end

filenirs = [ fpath  fpath2  fname ];

% LOAD DATA
mydat = load( [ fpath  fpath2  fname  ]  );
% OMIT UNUSED CHANS
%
% NOTES: 
% THERE IS A CORRECTION HERE FROM 20-20 TO 6-6
% ONLY 1-1, 2-2, 3-3, etc.
%
[nsamples, nchans0] = size( mydat );
xxnirs = zeros( nsamples, 19 );
for ichan = 1:nchans0,
    if 20==ichan, % we used 20-20 for chan 6
        xxnirs( : , 6 ) = mydat( : , ichan );
    else, % we used 1-1 for chan1, 2-2 for chan2 etc.
        xxnirs( : , ichan ) = mydat( : , ichan );
    end
end

%[nirschanlabels, xxnirsf ] = getpreprocess_nirs(filenirs,filterlofreq,filterhifreq, filternotch, Fsnirs);
% FILTER
xxnirsf = mybutterandiirnotchfilters( xxnirs, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
% NIRS CHAN LABELS
for ii = 1:19,
    nirschanlabels{ii} = get_chanlabel_from_srcdec( ii, ii );
    fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
end

%% GET NIRS EVENT FILE
fpath = 'F:\LabData\LabData\AhmetOmurtag1216\NIRS\';

if 1==icycle,
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_001 - eyesopen1\';
    fname = 'NIRS-2013-12-16_001'; 
elseif 2==icycle,
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_002 - eyesopen2\';
    fname = 'NIRS-2013-12-16_002';
elseif 3==icycle,
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_004 - eyesopen3\';
    fname = 'NIRS-2013-12-16_004';
elseif 4==icycle,
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_003 - eyesclose1\';
    fname = 'NIRS-2013-12-16_003'; 
elseif 5==icycle,
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_005 - eyesclose2\';
    fname = 'NIRS-2013-12-16_005';
elseif 6==icycle,
    fpath2 = 'Detector Readings\Ahmet Omurtag- 2013-12-16_006 - taskvft\';
    fname = 'NIRS-2013-12-16_006'; 
end
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
Tw = 60; % analysis window (sec)
Nwmax = 100000000; % max num windows to use (large number implies use max available)
nlags = ceil( Fsnirs.*(Tw-1) ); % max number of lags to compute (samples)

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

%%
figure(990); clf; ymaxlim = .1; 
pp1d = plot( Lags./Fsnirs, XCFmeansurr1d, 'g' ); hold on;
pp = plot( Lags./Fsnirs, XCFmean, 'r' ); hold on; 
axis([-maxlim  maxlim  -ymaxlim  ymaxlim]); grid on;
title(['eeg-nirs and eeg-eeglagged']);
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






