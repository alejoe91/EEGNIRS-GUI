clear;
%% PARAMS
% Input files

%    fpath = 'F:\LabData\LabData\AhmetOmurtag1216\EEG\';
%   fname = 'AhmetOmurtag-eyesclose1_20131216_193358.edf';
% fname = 'AhmetOmurtag-eyesopen1_20131216_191542.edf';
%  fname = 'AhmetOmurtag-closenotrigger_20131216_202346.edf';
% fname = 'AhmetOmurtag-eyesopen3_20131216_194123.edf';
%fname = 'AhmetOmurtag-taskvft_20131216_200735.edf';

fnames = {};

fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
fname = [fpath 'Omurtag327eyescloseonlyeeg_20140327_161725.edf'];
fnames = [fnames; fname];
fname = [fpath 'Omurtag327restingstateboth_20140327_162652.edf'];
fnames = [fnames; fname];
fname = [fpath 'Omurtag327restingstateeyesopenboth_20140327_163413.edf'];
fnames = [fnames; fname];



% PRINT FIGS
figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
flag_print_fig_conn = 0;
fnamefigconn = 'fig_eeg_conn_taskvft';
flag_print_fig_hist = 0;
fnamefighist = 'fig_eeg_hist_eyesopen1';

% ANALYSIS
Tw = 1; % analysis window (sec)
Nwmax = 100000000; % max num windows to use (large number implies use max available)
corrthreshold = .85; % correlation threshold used in binary adjacency matrix
corrlagthreshold = 2; % (samples) lag threshold used in binary adjacency matrix
nlags = 40; % max number of lags to compute (samples)
adjmatviewthreshold = 0.; % threshold for viewing connections in the final adj. matrix

% PREPROCESS
% data clip
dtbeg = 50; % sec
dtend = 5; % sec 
% Bandpass filter and notch (normal bandpass range: 0.5-70 Hz)
filterlofreq = .5; % Hz
filterhifreq = 70; % Hz
filternotch = 60; % 0: no notch; 
% Set this to n, to omit the top n PCAs from the signals
flag_pca_cleanup = 0;
% true: normalize each window data before using it for analysis. false:
% don't normalize
flag_normalize_datawindow = true;

% DISPLAY PARAMS
% figure number for Epoch averaged adjacency matrix
fignum_meanadjmat = 2042;

% visualize corr for each window
flag_viewperwindow = 0;

%  histogram of adjacency matrix avg over epoch
flag_plothist = true; % true: plot; false: don't plot
flag_plothistsemilogy = false; % true: plot semilog, false: plot linear

% skip this many when plotting the adjacency matrix avg over epochs
% set this to a large number in order to visualize only last one
nnw = 10000;

% max connector line thickness
maxconnectionlinewidth = 100;

% plot lagged correlation params
flag_plotlaggedcorr = false;
pausetime_plotlaggedcorr = .7;

% Calculate and plot time course of epoch-avg-adjacency-matrix
flag_plotadjacencymatrixevolution = true;

%% read and analyze all files in a loop
for ii = length(fnames),
    
    fname = char( fnames{ii} );
    fprintf('[%s]\n',fname);
    
%% READ DATA
numrecordstoread = 2000; % to read all records set numrecordstoread=[]
numchanmax = 26; %21;
[genheadr1, sigheadr1, totnumrecs1, xx, Fs, numchan] = readEDFfile2( fname, numrecordstoread, numchanmax );

%% CLIP BEGINNING AND END
nbeg = Fs.*dtbeg;
nend = Fs.*dtend;
xx = xx( nbeg:end-nend, : );
%% FILTER
xxf = mybutterandiirnotchfilters( xx, [filterlofreq filterhifreq filternotch 1], 6, Fs );
%xxf = mybandstopfilter( xxf, 6.25, 6, Fs );


%% initialize the names and positions of predetermined full set of 128 channels
[ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
subsetindxs = [];
fprintf('Channels Used:\n');
for ii = 1:numchan,
    chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
    if ifound1,
        subsetindxs = [subsetindxs; ii];           
        fprintf( '%d %s\n', ii, chan1);
    end
end
%% VIEW SIG
if 0
    nbeg = 1; %Fs.*1;
    nend = nbeg + Fs.*100;
    tt = [nbeg:nend]./Fs;
    for ii = 1:numchan,
        chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        if ifound1,
%            figure(110); clf; plot( tt, xxf(nbeg:nend, ii) ); title( [ num2str(ii) ':   ' chan1 ] ); 
            figure(1102); clf; pmtm( xxf(nbeg:nend, ii), 3.5, [], Fs); title( [ num2str(ii) ':   ' chan1 ] );
            pause(.5);
        else
            figure(112); clf; plot( tt, xxf(nbeg:nend, ii) ); title( [ num2str(ii) ':   ' chan1 ] ); 
            pause(.5);
        end
    end
end
%% PCA: for X where each row is an observation, and each column a variable,
%% x = a * E',  a = x * E
if flag_pca_cleanup
    myxxf = xxf(:, subsetindxs);
    Q = cov( myxxf ); [EE,EV] = eig(Q); ev = diag(EV); [ev,indx] = sort(ev,1,'descend');
    EE = EE(:,indx); figure(220); clf; plot( ev,'.-' ); title( 'PCA eigenvalues');
    aa = myxxf * EE;
end
%% VIEW PCA MODES
if 0
    nbeg = Fs.*1; nend = nbeg + Fs.*10;
    tt = [nbeg:nend]./Fs;
    for  imode = 1:19;
        figure(420); clf;
        plot( tt, aa(nbeg:nend, imode)  );
        title( num2str(imode), 'fontsize', 22 );
        pause(1)
    end
end
%% APPLY PCA CLEANUP
if flag_pca_cleanup>0,
    aa(:, 1:flag_pca_cleanup) = 0;
    myxxf = aa * EE';
    xxfnc(:, subsetindxs) = myxxf;
else
    xxfnc = xxf;
end
%% CLEAN UP (eg EYEBLINKS)
if 0
    xtmp = xxfnc( : , subsetindxs);
    xoc1 = xxfnc( : , 23 );
    xoc2 = xxfnc( : , 25 );
    [xtmpc] = remove_ocular_adapt_filter( xtmp, xoc1, xoc2, Fs ); % Required DSP toolbox
    xxfnc( : , subsetindxs ) = xtmpc;
end
%% VIEW SIG
if 0
    nbeg = 1; %Fs.*1;
    nend = nbeg + Fs.*10;
    tt = [nbeg:nend]./Fs;
    for ii = 1:numchan,
        chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        if ifound1,
            figure(110); clf; plot( tt, xxf(nbeg:nend, ii), tt, -0 + xxfnc(nbeg:nend, ii) ); title( [ num2str(ii) ':   ' chan1 ] ); 
            pause(.5);
        else
            figure(112); clf; plot( tt, xxf(nbeg:nend, ii), tt, -0 + xxfnc(nbeg:nend, ii) ); title( [ num2str(ii) ':   ' chan1 ] ); 
            pause(.5);
        end
    end
end
%% CORR OF FULL SERIES
if 0
    nSTDs = [];
    nLags = 1.*Fs;
    figure(330); clf;
    for ii = 1:numchan,
        for jj = ii:numchan,
            [XCF,Lags,Bounds] = crosscorr(  xxfn(:,ii), xxfn(:,jj),  nLags, nSTDs);
            plot( Lags, XCF ); title( [num2str(ii) '  ' num2str(jj)] ); axis([-nLags nLags -1 1]); grid on; pause(.25)
        end
    end
end
%% MAX OF THE LAGGED CORR WITHIN WINDOW
numsampW = Fs.*Tw;
nstds = [];

maxdim = min( Nwmax, floor( length(xxfnc)./numsampW ) );

% max correlation
corrmax = zeros( numchan, numchan, maxdim );
% lag at which corr is maximum
lagcorrmax = zeros( numchan, numchan, maxdim );

% first window
iw = 1;
ibeg = (iw-1).*numsampW + 1;
iend = ibeg + numsampW - 1;

fprintf( 'Calculating lagged correlation maxima for all windows all channel pairs\n' );

% loop through windows
while iw<=Nwmax & iend<=length(xxfnc),

    dat = xxfnc( ibeg:iend, : );
    if flag_normalize_datawindow, dat = mynormalize( dat ); end

    for ii = 1:numchan,
        chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        for jj = ii:numchan,
            chan2 = sigheadr1(jj).label; [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );

            if ifound1 & ifound2,
                %[ XCF, Lags, Bounds ] = crosscorr(  dat(:,ii), dat(:,jj),  nlags, nstds );
                [ XCF, Lags ] = xcorr(  dat(:,ii), dat(:,jj),  nlags, 'coeff' );
                
                [pks1,locs1] = findpeaks( XCF );
                [pks2,locs2] = findpeaks( -XCF );
                pks = [pks1; -pks2];
                locs = [locs1; locs2];
                [pksmaxdum, indxmax] = max( abs(pks) );
                pksmax = pks(indxmax);
                lagpksmax = Lags(locs(indxmax) );
                corrmax(ii, jj, iw) = pksmax;
                lagcorrmax(ii, jj, iw) = lagpksmax;
                
                if flag_plotlaggedcorr, % plot correlation vs. lag
                    figure(880); clf;
                    plot( Lags./Fs, XCF, 'k' ); hold on;
                    maxlim = max( Lags./Fs); axis([-maxlim maxlim -1 1]); grid on;
                    titlestr = [  'iw ' num2str(iw) ', chans ' chan1 ' ' chan2 ', ' num2str(ii) '  ' num2str(jj) ',  maxcorr ' num2str( corrmax(ii, jj, iw) ) ', lag maxcorr: ' num2str(lagcorrmax(ii, jj, iw)  ) ];
                    title( titlestr );
                    imid = (length(Lags)-1)./2 + 1;
                    plot( Lags(imid)./Fs, XCF(imid), 'bo' );
                    plot( Lags(locs)./Fs, pks, 'rs' );
                    plot( lagpksmax./Fs, pksmax, 'ko' );                    
                    plot( lagpksmax./Fs, pksmax, 'k.' );

                end
                
                
                if 0
                [cmxdum, indxmaxcorr] = max(  abs(XCF)  );
                cmx = XCF(indxmaxcorr);
                corrmax(ii, jj, iw) = cmx;
                lagcorrmax(ii, jj, iw) = Lags(indxmaxcorr);

                if flag_plotlaggedcorr, % plot correlation vs. lag
                    figure(880); clf;
                    plot( Lags./Fs, XCF, 'k' ); hold on;
                    maxlim = max( Lags./Fs); axis([-maxlim maxlim -1 1]); grid on;
                    titlestr = [  'iw ' num2str(iw) ', chans ' chan1 ' ' chan2 ', ' num2str(ii) '  ' num2str(jj) ',  maxcorr ' num2str( corrmax(ii, jj, iw) ) ', lag maxcorr: ' num2str(lagcorrmax(ii, jj, iw)  ) ];
                    title( titlestr );
                    imid = (length(Lags)-1)./2 + 1;
                    plot( Lags(imid)./Fs, XCF(imid), 'bo' );
                    plot( Lags(indxmaxcorr)./Fs, XCF(indxmaxcorr), 'rx' );
                    plot( Lags(indxmaxcorr)./Fs, cmx, 'gs' );
                    xlabel('Lag (s)'); ylabel('Correlation');
                    hold off;
                    pause(pausetime_plotlaggedcorr);
                end
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

%% VIEW FUNC NET FOR WINDOW WITH WEIGHTED CONNECTIONS
if flag_viewperwindow,
    gry = [.9 .9 .7];

    for iw = 1:Nw,
        fignum = 1020; figure(fignum);
        clf; hold on;
        title( {    [fname];    ['Corr within window'];...
            ['Window: ' num2str( Tw ) ', Thresh: ' num2str( corrthreshold ) ', ' num2str( corrlagthreshold ) ];...
            [ num2str( iw ) ' / ' num2str( Nw ) ]    }  ,'Interpreter', 'none' );

        for ii = 1:numchan,
            chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );

            for jj = ii:numchan,
                chan2 = sigheadr1(jj).label; [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );

                if ifound1 & ifound2 & corrlagthreshold<=abs(lagcorrmax(ii, jj, iw)) ...
                        & corrthreshold<=abs(corrmax(ii, jj, iw)),
                    mycorr = abs(corrmax(ii, jj, iw));
                    mylagcorr = lagcorrmax(ii, jj, iw);
                    fprintf( 'Window iw: %5d, chans: %6s %6s, corr %10.5f, corrlag %5d, Thresholds: corr  %10.5f, lag  %10.5f\n', ...
                        iw, chan1, chan2, mycorr, mylagcorr, corrthreshold, corrlagthreshold);
                    linesz = maxconnectionlinewidth.*mycorr;
                    plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)],'linewidth', linesz, 'color', gry );
                end
            end
        end
        iverbose = 0; labelclr = 'k'; rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr );
        pause( .001 );
    end
end
%% BUILD ADJACENCY MATRICES
% Binary adjacency matrix
AdjMatBin = zeros( numchan, numchan, maxdim );
% Adjacency matrix averaged over epochs (multiple windows, iw = 1,...,nw
AdjMatMean = zeros( numchan, numchan, maxdim );

fprintf( 'Calc adjacency matrices\n');
for iw = 1:Nw,
    for ii = 1:numchan,
        chan1 = sigheadr1(ii).label;[ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        for jj = ii:numchan,
            chan2 = sigheadr1(jj).label; [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
            if ifound1 & ifound2 & corrlagthreshold<=abs(lagcorrmax(ii, jj, iw)) ...
                    & corrthreshold<=abs(corrmax(ii, jj, iw)),
                AdjMatBin( ii, jj, iw ) = 1;
                fprintf( 'iw %5d / %5d, %6s %6s \n', iw, Nw, chan1, chan2 );
            end
        end
    end
    AdjMatMean( :, :, iw ) = mean( AdjMatBin, 3 );
end
%% HISTOGRAM OF EPOCH AVG CORR (ALL WINDOWS)
if flag_plothist,
    % flatten corr matrix
    [mm, nn, kk] = size( AdjMatMean );
    flatvector = reshape( AdjMatMean, mm.*nn.*kk, 1 );
    % remove zero elements
    k0 = find( 0==flatvector );
    flatvector(k0) = [];
    [pppp,xxxx] = hist( flatvector );
    figure( 7070); 
    if flag_plothistsemilogy, semilogy(xxxx,pppp,'.-'); 
    else, plot(xxxx,pppp,'.-');  end;
    xlabel( 'Epoch Averaged Correlation (All Windows)' ); ylabel( 'Frequency' ); title( 'Histogram (zeros removed)' );
end
%% HISTOGRAM OF EPOCH AVG CORR (FINAL WINDOW)
if flag_plothist,
    % flatten corr matrix
    [mm, nn, kk] = size( AdjMatMean );
    flatvector = reshape( AdjMatMean(:,:,end), mm.*nn, 1 );
    % remove zero elements
    k0 = find( 0==flatvector );
    flatvector(k0) = [];
    [pppp,xxxx] = hist( flatvector );
    figure( 7080); 
    if flag_plothistsemilogy, semilogy(xxxx,pppp,'.-'); 
    else, plot(xxxx,pppp,'.-');  end;
    xlabel( 'Epoch Averaged Correlation (Final Window)' ); ylabel( 'Frequency' ); title( 'Histogram (zeros removed)' );
end
%% VISUALIZE EPOCH AVG ADJ MATRIX
fh = figure(fignum_meanadjmat); set(fh,'color','w'); clf; hold on;
for iw = [ 1:nnw:Nw, Nw]
    clf;  hold on;
    title( {    fname;    ['Mean Adj. Matrix'];...
        [ 'Window size (s): ' num2str( Tw ) ', Thresh. corr: ' num2str( corrthreshold ) ', lag (samples): ' num2str( corrlagthreshold ) ', Corr view thresh: ' num2str( adjmatviewthreshold )];...
        [ num2str( iw ) ' / ' num2str( Nw ) ]    }  ,'Interpreter', 'none');
    for ii = 1:numchan,
        chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        for jj = ii:numchan,
            chan2 = sigheadr1(jj).label; [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
            if ifound1 & ifound2 & adjmatviewthreshold<AdjMatMean( ii, jj, iw ),
                fprintf( 'iw %5d / %5d, %6s %6s,  AdjMatMean %10.5f \n', iw, Nw, chan1, chan2,  AdjMatMean( ii, jj, iw ) );
                linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, iw ); % ./max(max(max(AdjMatMean)));
                gry = [.8 .9 .8];
                plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry );
            end
        end
    end
    iverbose = 0; labelclr = 'k'; rendermontage( fignum_meanadjmat, sigheadr1, labels, xchan, iverbose, labelclr );
    drawnow
end
%%
if flag_print_fig_conn,
    fignum = fignum_meanadjmat + 10000
    fh = figure(fignum ); set(fh,'color','w'); clf; hold on;
     for ii = 1:numchan,
        chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        for jj = ii:numchan,
            chan2 = sigheadr1(jj).label; [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
            if ifound1 & ifound2 & adjmatviewthreshold<AdjMatMean( ii, jj, iw ),
                fprintf( 'iw %5d / %5d, %6s %6s,  AdjMatMean %10.5f \n', iw, Nw, chan1, chan2,  AdjMatMean( ii, jj, iw ) );
                linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, iw ); % ./max(max(max(AdjMatMean)));
                gry = [.3 .3 .5]; %gry = [.8 .9 .8];
                plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry );
            end
        end
    end
    iverbose = 0; labelclr = 'k'; rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr );
    drawnow
    
    set( gcf, 'PaperPositionMode', 'auto')
    print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefigconn] );
end

%% ADJ MATRIX TIME EVOLUTION
if flag_plotadjacencymatrixevolution,
    
    [mm, nn, kk] = size( AdjMatMean );
    flatvectorfinal = reshape( AdjMatMean(:,:,Nw), mm.*nn, 1 );
    if 0==sum(flatvectorfinal), error( 'ZERO ADJACENCY MATRIX'); end
    ccchistory = []; iwhistory = [];
    for iw = [ 1:Nw ],
        flatvector = reshape( AdjMatMean(:,:,iw), mm.*nn, 1 );
        if 0==sum(flatvector), continue; end
        ccc = corr( flatvectorfinal, flatvector);%        fprintf('%d %f\n',iw,ccc);
        ccchistory = [ccchistory; ccc]; iwhistory = [iwhistory; iw];
    end

    figure(8080);clf;
    tccchistory =  iwhistory.*Tw;
    plot( tccchistory, ccchistory,'-'); axis([0 Tw*Nw 0 1]);
    xlabel('Time (s)'); ylabel('Correlation between current and final adjacency matrices');
 
   % save( 'ccchistory_eeg_vft.mat', 'tccchistory', 'ccchistory' )
       % save( 'ccchistory_eeg_eyesclose1.mat', 'tccchistory', 'ccchistory' )
    if flag_print_fig_hist,
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefighist] );
    end
end

%%
end % end read and analyze each files in a loop







