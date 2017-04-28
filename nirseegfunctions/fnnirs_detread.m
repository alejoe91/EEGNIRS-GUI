clear; fprintf('\n');
%{
Structure of this code

Assign input data filename
Read channel label corresponding to srcdet pair
Set param values
Load detector data
Convert to Hb
Preprocess
Calc correlation matrices
Calc adjacency matrices by thresholding

%}

for icycle = [301],
    %% PARAM
    
    if 9==icycle,
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\RestingState\';
        fname = 'NIRS-2014-03-27_002';
        badchans = {}; % example: badchans = {'FP1', 'T6' 'F4' 'FP1' 'F1' 'F8' 'P3' 'EKG1' };        
    elseif 10==icycle,
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\RestingStateeyesopen\';
        fname = 'NIRS-2014-03-27_003';
        badchans = {};
    elseif 11==icycle,
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\Taskeyesclose\';
        fname = 'NIRS-2014-03-27_004';
        badchans = {'CZ'};
    elseif 12==icycle,
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\Taskeyesopen\';
        fname = 'NIRS-2014-03-27_005';
        badchans = {};
    elseif 301==icycle,
        for ii = 1:20, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0416\NIRS\';
        fpath2 = 'Detector Readings\task\';
        fname = 'NIRS-2014-04-16_003';
        badchans = {};
    elseif 601==icycle,
        for ii = 1:20,
            chanlabel{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, chanlabel{ii} );
        end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0509\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-05-09_003';
        badchans = {'F8' 'T8' 'P8'};
    end
    % GET sample rate
    Fs = str2num( getfieldval( [fpath fpath2 fname '_config.txt'], 'SamplingRate', '=') );
    
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
    adjmatviewthreshold = 0.3;
    
    flag_normalize_datawindow = false;
    flag_plotadjacencymatrixevolution = true;
    %% DISPLAY PARAMS AND UTILS
    %figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
    figspath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
    
    YLEN = 450;
    subsetindxs = [1:19];
    idisplaytitle = 1;
    iprintfig = 0;
    iprintfig2 = 0;
    
    chanlabels0 = {};
    for ii=1:19,
        channame = get_chanlabel_from_srcdetpair( ii, ii );
        chanlabels0 = [chanlabels0; channame];
    end
    
    fignamestub = '_fnnirs';
    fignamestrn = [ '_bandp' num2str(filterlofreq) '-' num2str(filterhifreq) '_Tw' num2str(Tw) '_corrthresh' num2str( corrthreshold ) '_lagthresh' num2str( corrlagthreshold ) ];
    %% initialize the names and positions of predetermined full set of 128 channels
    [ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
    
    %% LOAD DATA
    mydat = load( [ fpath  fpath2  fname '.wl1' ]  );
    xx = fnnirs_reduce_date( mydat );
    
    % CORRECTION: chan 20 -> 17
    xx(:,17) = xx(:,20); xx(:,20) = [];
    xx_wl1 = xx;
    
    mydat = load( [ fpath  fpath2  fname '.wl2' ]  );
    xx = fnnirs_reduce_date( mydat );
    
    % CORRECTION: chan 20 -> 17
    xx(:,17) = xx(:,20); xx(:,20) = [];
    xx_wl2 = xx;
    %% CONVERT TO HEMOGLOBIN CONCENTRATION
    [xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
    
    %%
    for indx_oxy_deoxy = 1:2,
        
        if 1==indx_oxy_deoxy, % Oxy most chans GOOD
            fignamestub = '_fnnirs_HbO_';
            xxnirs = xxo;
            titlestrn = {...
                [fpath fpath2];...
                [fname fignamestub];...
                [ 'Window size (s): ' num2str( Tw ) ', Thresh. corr: ' num2str( corrthreshold ) ', lag (samples): ' num2str( corrlagthreshold ) ', Corr view thresh: ' num2str( adjmatviewthreshold )];...
                };
            
        else, % R deoxy most channels BAD
            fignamestub = '_fnnirs_HbR_';
            xxnirs = xxr;
            titlestrn = {...
                [fpath fpath2];...
                [fname fignamestub];...
                [ 'Window size (s): ' num2str( Tw ) ', Thresh. corr: ' num2str( corrthreshold ) ', lag (samples): ' num2str( corrlagthreshold ) ', Corr view thresh: ' num2str( adjmatviewthreshold )];...
                };
            
        end
        
        %% CLIP BEGINNING AND END
        nbeg = ceil(Fs.*dtbeg);
        nend = ceil(Fs.*dtend);
        xxnirs = xxnirs( nbeg:end-nend, : );
        %% FILTER
        xxnirsf = mybutterandiirnotchfilters( xxnirs, [filterlofreq filterhifreq filternotch 1], 6, Fs );
        %% VIEW DATA
        if 0
            for ii = 1:19,
                figure(330);clf;
                tt = [1:length(xxnirsf)]./Fs;
                plot( tt, xxnirsf(:, ii) );
                title( [num2str(ii)] , 'fontsize', 22);
                figure(440);clf; nw = 3.5; pmtm( xxnirs(:, ii),  nw, [], Fs);title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
                figure(450);clf;
                nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fs);
                title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
                figure(460);clf;
                nw = 3.5;
                [Pxx,F] = pmtm( xxnirsf(:, ii),  nw, [], Fs);
                logPxx = 10.*log10(Pxx);
                indxs = find( 0<F & F<.3 );
                plot(F,logPxx); axis([0 .3 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
                pause(.1);
            end
        end
        %% CLEANUP
        xxnirsfc = xxnirsf;
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
                            [pks1,locs1] = myfindpeaks( XCF );
                            [pks2,locs2] = myfindpeaks( -XCF );
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
        
        %% Individual figs
        if 1==iprintfig,
            
            %% PLOT ADJACENCY MATRICES
            fontsz=8;
            fh = figure(1000+icycle+indx_oxy_deoxy); set(fh,'color','w'); set(fh, 'position',[50 200 550 YLEN ]); clf; hold on;
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
            fignum = 4000 + icycle+indx_oxy_deoxy;
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
            fh = figure(8000+icycle+indx_oxy_deoxy); set(fh,'color','w'); set(fh, 'position',[950 200 350 YLEN]); clf; hold on;
            
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
            fh = figure(9000+icycle+indx_oxy_deoxy); set(fh,'color','w'); set(fh, 'position',[1300 200 350 YLEN]); clf; hold on;
            
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
        %% ALL FIGS IN ONE
        fignum = 10+icycle+indx_oxy_deoxy;
        fh = figure(fignum); set(fh,'color','w'); set(fh, 'position',[50 200 1600 YLEN ]); clf; hold on;
        % PLOT ADJACENCY MATRICES
        subplot('position',[0.03 0.12 .25 .74 ])
        
        fontsz=8;
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
        
        
        % HEAD CONNECTIVITY
        subplot('position',[0.34 0.12 .25 .74 ])
        
        for ii = 1:numchan,
            chan1 = get_chanlabel_from_srcdetpair( ii, ii ); [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
            for jj = ii:numchan,
                chan2 = get_chanlabel_from_srcdetpair( jj, jj ); [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
                if ifound1 & ifound2 & adjmatviewthreshold<AdjMatMean( ii, jj, end ),
                    linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, end ); % ./max(max(max(AdjMatMean)));
                    gry = [.3 .3 .5]; %gry = [.8 .9 .8];
                    plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry ); hold on;
                end
            end
        end
        gry = [.35 .35 .35];
        iverbose = 0; labelclr = 'k'; rendermontage2( fignum, labels, xchan, iverbose, labelclr );
        if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
        drawnow
        
        % CONVERGENCE OF CONNECTIVITY
        subplot('position',[0.7 0.12 .1 .74 ]), hold on
        
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
        
        tccchistory =  iwhistory.*Tw;
        plot( tccchistory, ccchistory,'-','linewidth',3,'color','k');
        axis([0 Tw*Nw 0 max(ccchistory)]);
        xlabel('Time (s)','fontsize',fontsz); ylabel('Convergence of Adjacency Matrix','fontsize',fontsz);
        set(gca,'fontsize',fontsz); box on;
        
        % HISTOGRAM OF CONNECTIVITY
        subplot('position',[0.86 0.12 .1 .74 ])
        
        % flatten corr matrix
        gry = [.4 .4 .4];
        [mm, nn, kk] = size( AdjMatMean );
        flatvector = reshape( AdjMatMean(:,:,end), mm.*nn, 1 );
        % remove zero elements
        k0 = find( 0==flatvector );
        flatvector(k0) = [];
        [pppp,xxxx] = hist( flatvector );
        
        % semilogy(xxxx,pppp,'.-');
        bar(xxxx,pppp, 'EdgeColor',gry,'FaceColor',gry); xlabel( 'Adjacency Weight','fontsize',fontsz ); ylabel( 'Frequency','fontsize',fontsz ); title( 'Histogram (zeros removed)' );
        axis([0 1 0 50])
        
        if iprintfig2,
            fnamefig =  ['fig_all' fignamestub fignamestrn '.tiff'];
            set( gcf, 'PaperPositionMode', 'auto')
            print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
        end
        
    end
    
end




