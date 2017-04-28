clear; fprintf('\n');

for icycle = [10],
    %% PARAM
    
    if 9==icycle,
        % MORE THAN HALF THE CHANNELS HAVE BAD DATA
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\RestingState\';
        fname = 'NIRS-2014-03-27_002_deoxyhb_T1to2919_C1to20.txt';
        badchans = {}; % example: badchans = {'FP1', 'T6' 'F4' 'FP1' 'F1' 'F8' 'P3' 'EKG1' };
    elseif 10==icycle,
        % CZ BAD
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\RestingState\';
        fname = 'NIRS-2014-03-27_002_oxyhb_T1to2919_C1to20.txt';
        badchans = {'CZ'};
    elseif 11==icycle,
        % TOO MANY BAD CHANNELS
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\RestingStateeyesopen\';
        fname = 'NIRS-2014-03-27_003_deoxyhb_T1to2238_C1to20.txt';
        badchans = {};
    elseif 12==icycle,
        % CZ PZ BAD
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\RestingStateeyesopen\';
        fname = 'NIRS-2014-03-27_003_oxyhb_T1to2238_C1to20.txt';
        badchans = {'CZ'};
    elseif 13==icycle,
        % TOO MANY CHANNELS HAVE BAD DATA
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\Taskeyesclose\';
        fname = 'NIRS-2014-03-27_004_deoxyhb_T1to2170_C1to20.txt';
        badchans = {};
    elseif 14==icycle,
        % CZ PZ BAD
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\Taskeyesclose\';
        fname = 'NIRS-2014-03-27_004_oxyhb_T1to2170_C1to20.txt';
        badchans = {};
    elseif 15==icycle,
        % TOO MANY CHANNELS HAVE BAD DATA
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\Taskeyesopen\';
        fname = 'NIRS-2014-03-27_005_deoxyhb_T1to2368_C1to20.txt';
        badchans = {};
    elseif 16==icycle,
        % CZ PZ P4 BAD
        Fs = 6.25;
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Hemoglobin Changes (NO FILTER)\Taskeyesopen\';
        fname = 'NIRS-2014-03-27_005_oxyhb_T1to2368_C1to20.txt';
        badchans = {};
        
    end
    %% PREPROCESS PARAM
    iflag_view_delayed_corr = false;
    % Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
    filterlofreq = .01; % Hz
    filterhifreq = 1.1; % Hz
    filternotch = 0; % 0: no notch;
    
    % PCA cleanup
    flag_pca_cleanup = 0;
    
    % data clip
    dtbeg = 10; % sec
    dtend = 10; % sec
    %% ANALYSIS PARAMS
    Tw = 4; % analysis window (sec)
    Nwmax = 100000000; % max num windows to use (large number implies use max available)
    corrthreshold = .9; % correlation threshold used in binary adjacency matrix
    corrlagthreshold = 0; % (samples) lag threshold used in binary adjacency matrix
    nlags = []; % max number of lags to compute (samples)
    numchan = 19;
    
    maxconnectionlinewidth = 20;
    nnw = 10000000;
    adjmatviewthreshold = 0.4;
    
    flag_normalize_datawindow = false;
    flag_plotadjacencymatrixevolution = true;
    %% DISPLAY PARAMS AND UTILS
    %figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
    figspath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
    
    YLEN = 450;
    subsetindxs = [1:19];
    idisplaytitle = 1;
    iprintfig = 0;
    
    chanlabels0 = {};
    for ii=1:19,
        channame = get_chanlabel_from_srcdetpair( ii, ii );
        chanlabels0 = [chanlabels0; channame];
    end
    
    titlestrn = {...
        fpath;...
        fname;...
        [ 'Window size (s): ' num2str( Tw ) ', Thresh. corr: ' num2str( corrthreshold ) ', lag (samples): ' num2str( corrlagthreshold ) ', Corr view thresh: ' num2str( adjmatviewthreshold )];...
        };
    
    fignamestub = '_fnnirs';
    fignamestrn = [ '_bandp' num2str(filterlofreq) '-' num2str(filterhifreq) '_Tw' num2str(Tw) '_corrthresh' num2str( corrthreshold ) '_lagthresh' num2str( corrlagthreshold ) ];
    
    %% LOAD DATA
    mydat = load( [ fpath  fpath2  fname  ]  );
    %% MAP NIRS SRC-DEC PAIR TO 10-20 CHAN LABEL
    %for ii = 1:19, chanlabel = get_chanlabel_from_srcdec( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
    %% initialize the names and positions of predetermined full set of 128 channels
    [ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
    %% OMIT UNUSED CHANS
    [nsamples, nchans0] = size( mydat );
    if 0
        %
        % NOTE:
        % FOR AhmetOmurtag1216
        % THERE IS A CORRECTION HERE FROM 20-20 TO 6-6
        %
        xxnirs = zeros( nsamples, 19 );
        for ichan = 1:nchans0,
            if 20==ichan, % we used 20-20 for chan 6
                xxnirs( : , 6 ) = mydat( : , ichan );
            else, % we used 1-1 for chan1, 2-2 for chan2 etc.
                xxnirs( : , ichan ) = mydat( : , ichan );
            end
        end
    end
    %
    % NOTE:
    % FOR AhmetOmurtag0327
    % THERE IS A CORRECTION HERE FROM 20-20 TO 17-17
    %
    xxnirs = zeros( nsamples, 19 );
    for ichan = 1:nchans0,
        if 20==ichan, % we used 20-20 for chan 17
            xxnirs( : , 17 ) = mydat( : , ichan );
        else, % we used 1-1 for chan1, 2-2 for chan2 etc.
            xxnirs( : , ichan ) = mydat( : , ichan );
        end
    end
    %% CLIP BEGINNING AND END
    nbeg = Fs.*dtbeg;    nend = Fs.*dtend;    xxnirs = xxnirs( nbeg:end-nend, : );
    
    %% VIEW DATA AND SPECTRUM
    if 0
        [nsamples, nchans] = size( xxnirs );
        for ichan = 1:nchans,
            series = xxnirs(:,ichan);
            tt = [1:length(xxnirsf)]./Fs;
            figure(110);clf;
            plot( tt, series ); title( [num2str(ichan) '  -  ' get_chanlabel_from_srcdetpair( ichan, ichan )  ] , 'fontsize', 22);
            figure(220);clf;
            nw = 3.5;
            pmtm( series,  nw, [], Fs); title( [num2str(ichan) '  -  ' get_chanlabel_from_srcdetpair( ichan, ichan )  ] , 'fontsize', 22);
            pause( 1 )
        end
    end
   
    %% FILTER
    xxnirsf = mybutterandiirnotchfilters( xxnirs, [filterlofreq filterhifreq filternotch 1], 6, Fs );
    %% VIEW DATA
    if 1
        for ii = 1:19,
            figure(330);clf;
            tt = [1:length(xxnirsf)]./Fs;
            plot( tt, xxnirsf(:, ii) );
            title( [num2str(ii)] , 'fontsize', 22);
            figure(440);clf; nw = 3.5; pmtm( xxnirs(:, ii),  nw, [], Fs);title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            figure(450);clf; nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fs);title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            pause(.1);
        end
    end
    
    %% PCA: for X where each row is an observation, and each column a variable,
    %% x = a * E',  a = x * E
    if flag_pca_cleanup
        Q = cov( xxnirsf ); [EE,EV] = eig(Q); ev = diag(EV); [ev,indx] = sort(ev,1,'descend');
        EE = EE(:,indx); figure(220); clf; plot( ev,'.-' ); title( 'PCA eigenvalues');
        aa = xxnirsf * EE;
        % VIEW PCA MODES
        if 0
            nbeg = Fs.*1; nend = nbeg + Fs.*240;
            tt = [nbeg:nend]./Fs;
            for  imode = 1:19;
                figure(420); clf;
                plot( tt, aa(nbeg:nend, imode) , tt, -.0001 + xxnirsf(nbeg:nend, imode) );
                title( num2str(imode), 'fontsize', 22 );
                figure(450);clf; nw = 3.5; pmtm( aa(nbeg:nend, imode),  nw, [], Fs);
                pause(.1)
            end
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
            chan1 = get_chanlabel_from_srcdetpair( ii, ii );
            indxbadchan1 = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
            if 0==indxbadchan1,
                for jj = ii+1:numchan,
                    chan2 = get_chanlabel_from_srcdetpair( jj, jj );
                    indxbadchan2 = findstrincellarray( badchans, upper(chan2) ); % return zero if chan1 is NOT a bad channel
                    if 0==indxbadchan2,
                        [ XCF, Lags ] = xcorr(  dat(:,ii), dat(:,jj),  nlags, 'coeff' );
                        [pks1,locs1] = findpeaks( XCF );
                        [pks2,locs2] = findpeaks( -XCF );
                        pks = [pks1; -pks2];
                        locs = [locs1; locs2];
                        [pksmaxdum, indxmax] = max( abs(pks) );
                        pksmax = pks(indxmax);
                        lagpksmax = Lags(locs(indxmax) );
                        if isempty(pksmax),
                            corrmax(ii, jj, iw) = 0;
                            lagcorrmax(ii, jj, iw) = 0;
                        else
                            corrmax(ii, jj, iw) = pksmax;
                            lagcorrmax(ii, jj, iw) = lagpksmax;
                        end
                        
                        if iflag_view_delayed_corr, % plot correlation vs. lag
                            figure(880); clf;
                            plot( Lags./Fs, XCF, 'k' ); hold on;
                            maxlim = max( Lags./Fs); axis([-maxlim maxlim -1 1]); grid on;
                            chan1 = get_chanlabel_from_srcdetpair( ii, ii );
                            chan2 = get_chanlabel_from_srcdetpair( jj, jj );
                            titlestr = [  'iw ' num2str(iw) ', chans ' chan1 ' ' chan2 ', ' num2str(ii) '  ' num2str(jj) ',  maxcorr ' num2str( corrmax(ii, jj, iw) ) ', lag maxcorr: ' num2str(lagcorrmax(ii, jj, iw)  ) ];
                            title( titlestr );
                            imid = (length(Lags)-1)./2 + 1;
                            plot( Lags(imid)./Fs, XCF(imid), 'bo' );
                            plot( Lags(locs)./Fs, pks, 'rs' );
                            plot( lagpksmax./Fs, pksmax, 'ko' );
                            plot( lagpksmax./Fs, pksmax, 'k.' );
                            pause(.1)
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
    %% BUILD ADJACENCY MATRICES
    % Binary adjacency matrix
    AdjMatBin = zeros( numchan, numchan, maxdim );
    % Adjacency matrix averaged over epochs (multiple windows, iw = 1,...,nw
    AdjMatMean = zeros( numchan, numchan, maxdim );
    
    fprintf( 'Calc adjacency matrices\n');
    for iw = 1:Nw,
        for ii = 1:numchan,
            chan1 = get_chanlabel_from_srcdetpair( ii, ii );
            [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
            for jj = ii+1:numchan,
                chan2 = get_chanlabel_from_srcdetpair( jj, jj );
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
    
    
    %% MOMA
    if 0 % "MOMA plot"
        if 4==icycle,
            [mm, nn, kk] = size( AdjMatMean );
            flatvector = reshape( AdjMatMean(:,:,Nw), mm.*nn, 1 );
            [flatvector,indx0] = sort(flatvector,'descend');
            figure(9020); clf; hold on;
            plot( smooth( log10(flatvector) ),'r');
        else
            flatvector = reshape( AdjMatMean(:,:,Nw), mm.*nn, 1 );
            flatvector = flatvector(indx0);
            flatvector = smooth( flatvector,20 );
            figure(9020);
            plot(  log10(flatvector),'g');
        end
    end
    
    %% PLOT ADJACENCY MATRICES
    fontsz=8;
    fh = figure(1000+icycle); set(fh,'color','w'); set(fh, 'position',[50 200 550 YLEN ]); clf; hold on;
    matrixplot0 = AdjMatMean(subsetindxs,subsetindxs,end);
    reindx = [ 1   9   3  17 2  16 10 4  18 6  15 11 5  19 7  14 12 8  13 ]; % reorder chans frontal, center, posterior
    chanlabels = chanlabels0 ( reindx, : );
    matrixplot = matrixplot0( reindx, reindx );
    matrixplot2 = matrixplot;
    for ii=1:length(matrixplot),
        for jj=ii+1:length(matrixplot),
            matrixplot2(jj,ii) = matrixplot(ii,jj);
        end
    end
    imagesc(  flipud(matrixplot2) ); %view([0,0,90]); colorbar
    mylen = length(chanlabels);
    for kk=1:mylen,
        text( kk-.3, mylen+1, upper( chanlabels(kk,:) ), 'fontsize', fontsz );
        text( -.4, mylen-kk+1, upper( chanlabels(kk,:) ), 'fontsize', fontsz );
        text( mylen+.8, mylen-kk+1, upper( chanlabels(kk,:) ), 'fontsize', fontsz );
    end
    axis off
    %colormapstrn = 'gray';
    colormap(1-gray);
    %        mymaxval = max(max(matrixplot));
    caxis([0, .5]);
    ch = colorbar('location','SouthOutside','fontsize',14);
    set(ch,'TickLength',[0 0])
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    cpos=get(ch,'Position');
    cpos(4)=cpos(4)/2; % Halve the thickness%% CALC TOP CONNECTED CHANS AND CHAN PAIRS
    cpos(3)=cpos(3).*.8; % Halve the thickness%% CALC TOP CONNECTED CHANS AND CHAN PAIRS
    cpos(2)=cpos(2) - 0.1; % Move it down
    cpos(1)=cpos(1) + 0.07; % Move it down
    set(ch,'Position',cpos);
    
    if iprintfig,
        fnamefig =  ['fig_connmatrix' fignamestub fignamestrn '.tiff'];
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
    end
    
    %% HEAD CONNECTIVITY
    fignum = 4000 + icycle;
    fh = figure(fignum ); set(fh,'color','w'); clf; hold on;set(fh, 'position',[600 200 350 YLEN]);
    
    for ii = 1:numchan,
        chan1 = get_chanlabel_from_srcdetpair( ii, ii ); [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        for jj = ii:numchan,
            chan2 = get_chanlabel_from_srcdetpair( jj, jj ); [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
            if ifound1 & ifound2 & adjmatviewthreshold<AdjMatMean( ii, jj, end ),
                linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, end ); % ./max(max(max(AdjMatMean)));
                gry = [.3 .3 .5]; %gry = [.8 .9 .8];
                plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry );
            end
        end
    end
    gry = [.35 .35 .35];
    iverbose = 0; labelclr = 'k'; rendermontage2( fignum, labels, xchan, iverbose, labelclr );
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    drawnow
    
    if iprintfig,
        fnamefig =  ['fig_connhead' fignamestub fignamestrn '.tiff'];
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
    end
    
    %% CONVERGENCE OF CONNECTIVITY
    [mm, nn, kk] = size( AdjMatMean );
    ccchistory = []; iwhistory = [];
    matfinal = AdjMatMean(:,:,end);
    flatvectorfinal = reshape( matfinal, mm.*nn, 1 );
    
    for iw = [ 1:Nw ],
        mat1 = AdjMatMean(:,:,iw);
        flatvector1 = reshape( mat1, mm.*nn, 1 );
        ccc = corr( flatvector1, flatvectorfinal );
        ccchistory = [ccchistory; ccc]; iwhistory = [iwhistory; iw];
    end
    
    fontsz=14;
    fh = figure(8000+icycle); set(fh,'color','w'); set(fh, 'position',[950 200 350 YLEN]); clf; hold on;
    
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    tccchistory =  iwhistory.*Tw;
    plot( tccchistory, ccchistory,'-','linewidth',3,'color','k');
    axis([0 Tw*Nw 0 max(ccchistory)]);
    xlabel('Time (s)','fontsize',fontsz); ylabel('Convergence of Adjacency Matrix','fontsize',fontsz);
    set(gca,'fontsize',fontsz); box on;
    
    
    if iprintfig,
        fnamefig =  ['fig_connconverge' fignamestub fignamestrn '.tiff'];        
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
    end
    
    %% HISTOGRAM OF CONNECTIVITY
    
    % flatten corr matrix
    gry = [.4 .4 .4];
    [mm, nn, kk] = size( AdjMatMean );
    flatvector = reshape( AdjMatMean(:,:,end), mm.*nn, 1 );
    % remove zero elements
    k0 = find( 0==flatvector );
    flatvector(k0) = [];
    [pppp,xxxx] = hist( flatvector );
    fh = figure(9000+icycle); set(fh,'color','w'); set(fh, 'position',[1300 200 350 YLEN]); clf; hold on;
    
    % semilogy(xxxx,pppp,'.-');
    bar(xxxx,pppp, 'EdgeColor',gry,'FaceColor',gry); xlabel( 'Adjacency Weight','fontsize',fontsz ); ylabel( 'Frequency','fontsize',fontsz ); title( 'Histogram (zeros removed)' );
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    axis([0 1 0 50])
    
    if iprintfig,
        fnamefig =  ['fig_connhistogram' fignamestub fignamestrn '.tiff'];
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
    end
    
    
    
end




