clear;
%% PARAMS
% Input files
fpath = 'C:\downloads\edffile\';
%fname = 'Ranjitsingh_Irwin_20121029_155839.edf';
%fname = 'pp.edf';
%fname = 'f62.edf';
% fname = 'AhmetOmurtag-eyesclose1_20131216_193358.edf';
%fname = 'AhmetOmurtag-eyesopen1_20131216_191542.edf';
%  fname = 'AhmetOmurtag-closenotrigger_20131216_202346.edf';
% fname = 'AhmetOmurtag-eyesopen3_20131216_194123.edf';
 fname = 'AhmetOmurtag-taskvft_20131216_200735.edf';

Tw = 1; % analysis window (sec)
Nwmax = 100000000; % max num windows to use
corrthreshold = .7; % correlation threshold used in adjacency matrix
corrlagthreshold = 2; % (samples) lag threshold used in adjacency matrix
nlags = 40; % max number of lags to compute (samples)
adjmatviewthreshold = 0.;
% Bandpass filter and notch (normal clinical bandpass range is 0.5-70 Hz)
filterlofreq = .5; % Hz
filterhifreq = 70; % Hz
filternotch = 60; % 0: no notch; 

% figure number for Epoch averaged adjacency matrix
fignum_meanadjmat = 2040;

% visualize corr for each window
flag_viewperwindow = 0;

%  histogram of adjacency matrix avg over epoch
flag_plothist = true; % true: plot; false: don't plot
flag_plothistsemilogy = false; % true: plot semilog, false: plot linear

% skip this many when plotting the adjacency matrix avg over epochs
% set this to a large number in order to visualize only last one
nnw = 10000;

% true: normalize each window data before using it for analysis. false: don't normalize
flag_normalize_datawindow = true;

% max connector line thickness
maxconnectionlinewidth = 12;

% plot lagged correlation params
flag_plotlaggedcorr = false;
pausetime_plotlaggedcorr = .7;

% Calculate and plot time course of epoch-avg-adjacency-matrix
flag_plotadjacencymatrixevolution = true;

%% READ DATA
numrecordstoread = 2000; % to read all records set numrecordstoread=[]
numchanmax = 21;
file1 = [fpath fname];
[genheadr1, sigheadr1, totnumrecs1, xx, Fs, numchan] = readEDFfile2( file1, numrecordstoread, numchanmax );

%% CLIP
dtbeg = 5; % sec
dtend = 5; % end
nbeg = Fs.*dtbeg;
nend = Fs.*dtend;
xx = xx( nbeg:end-nend, : );

%% FILTER
xxf = mybutterandiirnotchfilters( xx, [filterlofreq filterhifreq filternotch 1], 6, Fs );

%% NORMALIZE
% xxfn = mynormalize( xxf );
xxfn = xxf;

%% VIEW SIG
if 0
    nbeg = 1; %Fs.*1;
    nend = nbeg + Fs.*20;
    tt = [nbeg:nend]./Fs;
    for ichan = 1:19,
        figure(110); clf; plot( tt, xxf(nbeg:nend, ichan) ); title( num2str(ichan) ); pause(1);
    end
end
%% PCA: for X where each row is an observation, and each column a variable,
%% COV(X) is the covariance matrix.
if 0
    Q = cov( xxfn ); [EE,EV] = eig(Q); ev = diag(EV); [ev,indx] = sort(ev,1,'descend');
    EE = EE(:,indx); figure(220); clf; plot( ev,'.-' );
    aa = xxfn * EE;
    %% VIEW PCA MODES
    nbeg = Fs.*19; nend = nbeg + Fs.*20;
    tt = [nbeg:nend]./Fs;
    imode = 5;
    figure(420); clf; plot( tt, aa(nbeg:nend, imode) );
end
%% ICA
if 0,
    % Row j of W is the jth eigenvector.
    % s=Wx and x=As
    [icasig, A, W] = fastica ( xxfn' );
    icasig = icasig';
    %% VIEW ICA
    nbeg = Fs.*19; nend = nbeg + Fs.*20;
    tt = [nbeg:nend]./Fs;
    figure(550); clf;
    for imode = 1:numchan,
        plot( tt, icasig(nbeg:nend, imode) ); title( num2str(imode) );
        pause ( 1 );
    end
    %%
    iaa = xxfn * W';
    %%
    figure(560); clf;
    imode = 2,
    plot( tt, iaa(nbeg:nend, imode) ); title( num2str(imode) );
    %% REMOVE EYEBLINK
    imoderemove = [2, 10];
    Wc = W;
    Wc(imoderemove, :) = 0;
    xxfnc = ( xxfn * Wc' ) * A';
    %% VIEW CLEAN SIGNAL
    ichan = 2;
    figure(660); clf; hold on;
    plot( tt, xxfn(nbeg:nend, ichan), 'b' );
    plot( tt, xxfnc(nbeg:nend, ichan), 'r' );
end
%%
xxfnc = xxfn;
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
%% initialize the names and positions of predetermined full set of 128 channels
[ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
fprintf('Channels Used:\n');
for ii = 1:numchan,
    chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
    if ifound1,
        fprintf( '%d %s\n',ii,chan1);
    end
end
%% MAX OF THE CORR WITHIN PER WINDOW
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
                [ XCF, Lags, Bounds ] = crosscorr(  dat(:,ii), dat(:,jj),  nlags, nstds );
                [cmx, indxmaxcorr] = max(XCF);
                corrmax(ii, jj, iw) = cmx;
                lagcorrmax(ii, jj, iw) = Lags(indxmaxcorr);

                if flag_plotlaggedcorr, % plot correlation vs. lag
                    figure(880); clf;
                    plot( Lags, XCF, 'k' ); hold on;
                    maxlim = max(Lags); axis([-maxlim maxlim -1 1]); grid on;
                    titlestr = [  'iw ' num2str(iw) ', chans ' chan1 ' ' chan2 ', ' num2str(ii) '  ' num2str(jj) ',  maxcorr ' num2str( corrmax(ii, jj, iw) ) ', lag maxcorr: ' num2str(lagcorrmax(ii, jj, iw)  ) ];
                    title( titlestr );
                    imid = (length(Lags)-1)./2 + 1;
                    plot( Lags(imid), XCF(imid), 'bo' );
                    plot( Lags(indxmaxcorr), XCF(indxmaxcorr), 'rx' );
                    plot( Lags(indxmaxcorr), cmx, 'gs' );
                    hold off;
                    pause(pausetime_plotlaggedcorr);
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

                if ifound1 & ifound2 & corrlagthreshold<abs(lagcorrmax(ii, jj, iw)) ...
                        & corrthreshold<abs(corrmax(ii, jj, iw)),
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
            if ifound1 & ifound2 & corrlagthreshold<abs(lagcorrmax(ii, jj, iw)) ...
                    & corrthreshold<abs(corrmax(ii, jj, iw)),
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
                linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, iw )./max(max(max(AdjMatMean)));
                gry = [.8 .9 .8];
                plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry );
            end
        end
    end
    iverbose = 0; labelclr = 'k'; rendermontage( fignum_meanadjmat, sigheadr1, labels, xchan, iverbose, labelclr );
    drawnow
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
    plot(iwhistory.*Tw, ccchistory,'-'); axis([0 Tw*Nw 0 1]);
    xlabel('Time (s)'); ylabel('Correlation between current and final adjacency matrices');
end

%%







