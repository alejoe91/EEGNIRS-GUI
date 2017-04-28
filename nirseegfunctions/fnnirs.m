clear;
%% PARAM

%fpath = 'C:\downloads\nirseegdata\LabData\LabData\AhmetOmurtag1216\NIRS\';
fpath = 'F:\LabData\LabData\AhmetOmurtag1216\NIRS\';

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

% PRINT FIG NAMES
figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
flag_print_fig_conn = 0;
fnamefigconn = 'fig_nirs_conn_eyesopen1';
flag_print_fig_hist = 0;
fnamefighist = 'fig_nirs_hist_eyesopen1';

% ANALYSIS
Tw = 4; % analysis window (sec)
Nwmax = 100000000; % max num windows to use (large number implies use max available)
corrthreshold = .8; % correlation threshold used in binary adjacency matrix
corrlagthreshold = 0; % (samples) lag threshold used in binary adjacency matrix
nlags = []; % max number of lags to compute (samples)
numchan = 19;

maxconnectionlinewidth = 10;
fignum_meanadjmat = 6060;
nnw = 10000000;
adjmatviewthreshold = 0.;

flag_normalize_datawindow = false;

flag_plotadjacencymatrixevolution = true;

flag_plothist = true;
flag_plothistsemilogy = false;

% Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
filterlofreq = .01; % Hz
filterhifreq = .8; % Hz
filternotch = 0; % 0: no notch;

% PCA cleanup 
flag_pca_cleanup = 1;

% Sample rate
Fs = 3.125;

%% LOAD DATA
mydat = load( [ fpath  fpath2  fname  '.'  fnamext]  );
%% MAP NIRS SRC-DEC PAIR TO 10-20 CHAN LABEL
for ii = 1:19,
    chanlabel = get_chanlabel_from_srcdec( ii, ii );
    fprintf( '%d %d  %s\n', ii, ii, chanlabel );
end
%% initialize the names and positions of predetermined full set of 128 channels
[ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
%% VIEW DATA AND SPECTRUM
if 0
    [nsamples, nchans] = size( mydat );
    for isrc = 1:20,
        for idet = 1:20,
            if isrc~=idet, continue; end
            ichan = (isrc - 1).*20 + idet;
            series = mydat(:,ichan);
            tt = [1:nsamples]./Fs;
            figure(110);clf;
            plot( tt, series ); title( [num2str(isrc) ' ' num2str(idet) ] , 'fontsize', 22);
            figure(220);clf;
            nw = 3.5;
            pmtm( series,  nw, [], Fs);
            pause( .1 )
            %    if isrc==idet, pause; end
        end
    end
end
%% REDUCE DATA - OMIT UNUSED CHANS
%
% NOTES: 
% THERE IS A CORRECTION HERE FROM 20-20 TO 6-6
% ONLY 1-1, 2-2, 3-3, etc.
%
[nsamples, nchans0] = size( mydat );
xxnirs = zeros( nsamples, 19 );
for isrc = 1:20,
    for idet = 1:20,
        if isrc==idet, % we used 1-1, 2-2, etc.
            ichan = (isrc - 1).*20 + idet;
            if 20==isrc, % we used 20-20 for chan 6
                xxnirs( : , 6 ) = mydat( : , ichan );
            else, % we used 1-1 for chan1, 2-2 for chan2 etc.
                xxnirs( : , isrc ) = mydat( : , ichan );
            end
        end
    end
end
%% FILTER
xxnirsf = mybutterandiirnotchfilters( xxnirs, [filterlofreq filterhifreq filternotch 1], 6, Fs );
%% VIEW DATA
if 0
    for ii = 1:19,
        figure(330);clf;
        tt = [1:nsamples]./Fs;
        plot( tt, xxnirsf(:, ii) );
        title( [num2str(ii)] , 'fontsize', 22);
        figure(440);clf; nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fs);
        figure(450);clf; nw = 3.5; pmtm( xxnirs(:, ii),  nw, [], Fs);
        pause(1);
    end
end
%% PCA: for X where each row is an observation, and each column a variable,
%% x = a * E',  a = x * E
if flag_pca_cleanup
    Q = cov( xxnirsf ); [EE,EV] = eig(Q); ev = diag(EV); [ev,indx] = sort(ev,1,'descend');
    EE = EE(:,indx); figure(220); clf; plot( ev,'.-' ); title( 'PCA eigenvalues');
    aa = xxnirsf * EE;
end
%% VIEW PCA MODES
if 0
    nbeg = Fs.*1; nend = nbeg + Fs.*240;
    tt = [nbeg:nend]./Fs;
    for  imode = 1:19;
        figure(420); clf;
        plot( tt, aa(nbeg:nend, imode) , tt, -.2 + xxnirsf(nbeg:nend, imode) );
        title( num2str(imode), 'fontsize', 22 );
        pause(.1)
    end
end
%% APPLY PCA CLEANUP
if flag_pca_cleanup>0,
    if flag_pca_cleanup>0,
        aa(:, 1) = 0;
    elseif flag_pca_cleanup>1,
        aa(:, 2) = 0;
    elseif flag_pca_cleanup>2,
        aa(:, 3) = 0;
    elseif flag_pca_cleanup>3,
        aa(:, 4) = 0;
    end
    xxnirsfc = aa * EE';
else
    xxnirsfc = xxnirsf;
end
%% ICA
if 0
    % Row j of W is the jth eigenvector.
    % s=Wx and x=As
    [icasig, A, W] = fastica ( xxnirsf' );
    icasig = icasig';
    whos icasig
    %% VIEW ICA
    nbeg = floor( Fs.*1 ) ; nend = floor( nbeg + Fs.*240 );
    
    for imode = 1:16,
        tt = [nbeg:nend]./Fs; % tt = [1:length(icasig)]./Fs;
        figure(550); clf; plot(tt, icasig(nbeg:nend, imode) , tt, -.2 + xxnirsf(nbeg:nend, imode)  );
        title( num2str(imode), 'fontsize', 22 );
        figure(440);clf; nw = 3.5; pmtm( icasig(:, imode),  nw, [], Fs);
        pause(1)
    end
    
    icasig(:, 9) = 0;
    xxnirsfc = A * icasig;
end
%% VIEW DATA
if 0
    for ii = 1:19,
        figure(330);clf;
        tt = [1:nsamples]./Fs;
        plot( tt, xxnirsfc(:, ii),  tt, -.04 + xxnirsf(:, ii) );
        title( [num2str(ii)] , 'fontsize', 22);
        figure(440);clf; nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fs);
        figure(450);clf; nw = 3.5; pmtm( xxnirsfc(:, ii),  nw, [], Fs);
        pause(2);
    end
end

%% MAX OF THE LAGGED CORR WITHIN WINDOW

numsampW = floor( Fs.*Tw );
nstds = [];

maxdim = min( Nwmax, floor( length(xxnirsfc)./numsampW ) );

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
while iw<=Nwmax & iend<=length(xxnirsfc),

    dat = xxnirsfc( ibeg:iend, : );
    if flag_normalize_datawindow, dat = mynormalize( dat ); end

    for ii = 1:numchan,
        for jj = ii:numchan,
            [ XCF, Lags, Bounds ] = crosscorr(  dat(:,ii), dat(:,jj),  nlags, nstds );
            [cmx, indxmaxcorr] = max(XCF);
            corrmax(ii, jj, iw) = cmx;
            lagcorrmax(ii, jj, iw) = Lags(indxmaxcorr);

            if false, % plot correlation vs. lag
                figure(880); clf;
                plot( Lags./Fs, XCF, 'k' ); hold on;
                maxlim = max(Lags./Fs); axis([-maxlim maxlim -1 1]); grid on;
                titlestr = [  'iw ' num2str(iw) ', chans ' num2str(ii) '  ' num2str(jj) ',  maxcorr ' num2str( corrmax(ii, jj, iw) ) ', lag maxcorr: ' num2str(lagcorrmax(ii, jj, iw)  ) ];
                title( titlestr );
                imid = (length(Lags)-1)./2 + 1;
                plot( Lags(imid)./Fs, XCF(imid), 'bo' );
                plot( Lags(indxmaxcorr)./Fs, XCF(indxmaxcorr), 'rx' );
                plot( Lags(indxmaxcorr)./Fs, cmx, 'gs' );
                    xlabel('Lag (s)'); ylabel('Correlation');
                hold off;
                pause(.01);
            end
        end
    end    
    iw = iw + 1;
    % bounds of data for next window
    ibeg = (iw-1).*numsampW + 1;
    iend = ibeg + numsampW - 1;
end
Nw = iw - 1; % Total num of windows
%% BUILD ADJACENCY MATRICES
% Binary adjacency matrix
AdjMatBin = zeros( numchan, numchan, maxdim );
% Adjacency matrix averaged over epochs (multiple windows, iw = 1,...,nw
AdjMatMean = zeros( numchan, numchan, maxdim );

fprintf( 'Calc adjacency matrices\n');
for iw = 1:Nw,
    for ii = 1:numchan,
        chan1 = get_chanlabel_from_srcdec( ii, ii );
        [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );        
        for jj = ii:numchan,
            chan2 = get_chanlabel_from_srcdec( jj, jj );
            [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
            if ifound1 & ifound2 & corrlagthreshold<=abs(lagcorrmax(ii, jj, iw)) ...
                    & corrthreshold<=abs(corrmax(ii, jj, iw)),
                AdjMatBin( ii, jj, iw ) = 1;
                fprintf( 'iw %5d / %5d, %6s %6s \n', iw, Nw, chan1, chan2 );
            end
        end
    end
    AdjMatMean( :, :, iw ) = mean( AdjMatBin, 3 );
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
        chan1 = get_chanlabel_from_srcdec( ii, ii );
        [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        for jj = ii:numchan,
            chan2 = get_chanlabel_from_srcdec( jj, jj );
            [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
            if ifound1 & ifound2 & adjmatviewthreshold<AdjMatMean( ii, jj, iw ),
                fprintf( 'iw %5d / %5d, %6s %6s,  AdjMatMean %10.5f \n', iw, Nw, chan1, chan2,  AdjMatMean( ii, jj, iw ) );
                linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, iw )./max(max(max(AdjMatMean)));
                gry = [.8 .9 .8];
                plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry );
            end
        end
    end
    iverbose = 0; labelclr = 'k'; rendermontage2( fignum_meanadjmat, labels, xchan, iverbose, labelclr );
    drawnow
end

if flag_print_fig_conn,
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
    plot( iwhistory.*Tw, ccchistory,'-' ); axis([0 Tw*Nw 0 1]);
    xlabel('Time (s)'); ylabel('Correlation between current and final adjacency matrices');
    if flag_print_fig_hist,
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefighist] );
    end
end
%%










