clear; format compact; fprintf('\n');
%% DATA FILE PARAMS

%  path
parentpath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\corrmat\';

% Fill in the arrays for files
% Processes files that match condition (str1 & str2 & (~fexclude_str))
finclude_str1 = 'Omurtag327';
finclude_str2 = '';
fexclude_str = '';

%% MAP PARAMS
corrthreshold = .85; % correlation threshold used in binary adjacency matrix
corrlagthreshold = 2; % (samples) lag threshold used in binary adjacency matrix
nlags = 40; %

% max connector line thickness
maxconnectionlinewidth = 100;
% threshold for plotting a correlation on head;
adjmatviewthreshold = 0.01;
%% read correlation matrices and analyze

allFiles = dir( parentpath );
allNames = { allFiles.name };

for ifile = 3:length(allNames),
    fprintf('[%s]\n',allNames{ifile});
    fname = char( allNames{ifile} );
    
    if ( isempty(finclude_str1) | strfind(fname, finclude_str1) ) ...
            & ( isempty(finclude_str2) | strfind(fname, finclude_str2) ) ...
            & ( isempty(fexclude_str) | isempty( strfind(fname, fexclude_str) ) )
        
        fullfname = [ parentpath  fname ];
        loadstrng = ['load ' fullfname];
        fprintf('   >>> [%s]\n', loadstrng);
        eval( loadstrng );
        %% BUILD ADJACENCY MATRICES
        fprintf( 'Calc adjacency matrices\n');
        [ numchan, ndummy, maxdim ] = size(corrmax);
        % initialize Binary adjacency matrix
        AdjMatBin = zeros( numchan, numchan, maxdim );
        % initialize Adjacency matrix averaged over epochs (multiple windows, iw = 1,...,nw
        AdjMatMean = zeros( numchan, numchan, maxdim );
        
        for iw = 1:Nw,
            for ii = 1:numchan,
                if isempty( intersect(subsetindxs, ii) ), ifound1=0; else, ifound1=1; end;
                for jj = ii:numchan,
                    if isempty( intersect(subsetindxs, jj) ), ifound2=0; else, ifound2=1; end;
                    if ifound1 & ifound2 & corrlagthreshold<=abs(lagcorrmax(ii, jj, iw)) ...
                            & corrthreshold<=abs(corrmax(ii, jj, iw)),
                        AdjMatBin( ii, jj, iw ) = 1;
                        chan1 = sigheadr1(ii).label; chan2 = sigheadr1(jj).label;
                        fprintf( 'iw %5d / %5d, %6s %6s \n', iw, Nw, chan1, chan2 );
                    end
                end
            end
            AdjMatMean( :, :, iw ) = mean( AdjMatBin, 3 );
        end
        %% PLOT ADJACENCY MATRICES
        fontsz=10;
        fh = figure(1000+ifile); set(fh,'color','w'); set(fh, 'position',[100 200 700 650 ]); clf; hold on;
        chanlabels0 = char(sigheadr1(  subsetindxs ).label);
        matrixplot0 = AdjMatMean(subsetindxs,subsetindxs,end);
        reindx = [ 7 8 3 4 9 10 17 1 2 11 12 18 5 6 13 14 19 15 16]; % reorder chans
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
        colormapstrn = 'default';
        colormap(colormapstrn);
        %        mymaxval = max(max(matrixplot));
        caxis([0, .1]);
        ch = colorbar('location','SouthOutside','fontsize',14);
        set(ch,'TickLength',[0 0])
        title( fname , 'Interpreter','none','color','b');
        fadsfdadf=789879;
        %% CALC TOP CONNECTED CHANS AND CHAN PAIRS
        %% HEAD CONNECTIVITY        
        fignum = 4000 + ifile;
        fh = figure(fignum ); set(fh,'color','w'); clf; hold on;set(fh, 'position',[700 200 550 650 ]);
        for ii = 1:numchan,
            chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
            for jj = ii:numchan,
                chan2 = sigheadr1(jj).label; [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
                if ifound1 & ifound2 & adjmatviewthreshold<AdjMatMean( ii, jj, end ),
                    linesz = maxconnectionlinewidth.*AdjMatMean( ii, jj, end ); % ./max(max(max(AdjMatMean)));
                    gry = [.3 .3 .5]; %gry = [.8 .9 .8];
                    plot( [xxx1(1) xxx2(1)], [xxx1(2) xxx2(2)], 'linewidth', linesz, 'color', gry );
                end
            end
        end
        iverbose = 0; labelclr = 'k'; rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr );        
        title( fname , 'Interpreter','none','color','b');
        drawnow
        
        % set( gcf, 'PaperPositionMode', 'auto')
        % print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefigconn] );
        
        %% HISTOGRAM OF CONNECTIVITY
        
        %% CONVERGENCE OF CONNECTIVITY
        
        
        
        fads=8;
        
    end
    
end

%% read and analyze
return;
