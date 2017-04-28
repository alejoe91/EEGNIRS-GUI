clear; format compact; fprintf('\n');
%% DATA FILE PARAMS
%  path
parentpath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\corrmat\';

% Processes files that match condition (finclude_str1 & finclude_str2 & (~fexclude_str1 & fexclude_str2))
finclude_str1 = 'Omurtag416';
finclude_str2 = 'task';
finclude_str3 = '0.5-70';
fexclude_str1 = '';
fexclude_str2 = '';

%% MAP PARAMS
corrthreshold = .75; % correlation threshold used in binary adjacency matrix
corrlagthreshold = 2; % (samples) lag threshold used in binary adjacency matrix
nlags = 40; %
% max connector line thickness
maxconnectionlinewidth = 100;
% threshold for plotting a correlation on head;
adjmatviewthreshold = 0.01;
%% FIGURE PARAM
figspath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
iprintfig = 0; % 1: print & save individual figs
iprintfig2 = 0; % 1: save 1x4 figure to disc
YLEN = 450;
substrng = ['fneeg_thresh' num2str(corrthreshold) '_lag' num2str(corrlagthreshold) ]
titlestrn = {...
    [ substrng '_' finclude_str1 '+' finclude_str2 '+' finclude_str3 '-' fexclude_str1 '-' fexclude_str2];...
    ['viewthresh' num2str( adjmatviewthreshold) ];...
    };
idisplaytitle = 1;
%% read correlation matrices and analyze

allFiles = dir( parentpath );
allNames = { allFiles.name };

count = 1;
for ifile = 3:length(allNames),
    fprintf('[%s]\n',allNames{ifile});
    fname = char( allNames{ifile} );
    
    if ( isempty(finclude_str1) | strfind(fname, finclude_str1) ) ...
            & ( isempty(finclude_str2) | strfind(fname, finclude_str2) ) ...
            & ( isempty(finclude_str3) | strfind(fname, finclude_str3) ) ...
            & ( isempty(fexclude_str1) | isempty( strfind(fname, fexclude_str1) ) ) ...
            & ( isempty(fexclude_str2) | isempty( strfind(fname, fexclude_str2) ) )
        
        fullfname = [ parentpath  fname ];
        loadstrng = ['load ' fullfname];
        fprintf('   >>> [%s]\n', loadstrng);
        eval( loadstrng );
        
        %% BUILD ADJACENCY MATRICES
        fprintf( 'Calc adjacency matrices\n');
        [ numchan, ndummy, maxdim ] = size(corrmax);
        % initialize Binary adjacency matrix
        AdjMatBin = zeros( numchan, numchan, maxdim );
        
        Nw = length(corrmax);
        
        for iw = 1:Nw,
            for ii = 1:numchan,
                if isempty( intersect(subsetindxs, ii) ), ifound1=0; else, ifound1=1; end;
                for jj = ii:numchan,
                    if isempty( intersect(subsetindxs, jj) ), ifound2=0; else, ifound2=1; end;
                    if ifound1 & ifound2 & corrlagthreshold<=abs(lagcorrmax(ii, jj, iw)) ...
                            & corrthreshold<=abs(corrmax(ii, jj, iw)),
                        AdjMatBin( ii, jj, iw ) = 1;
                        chan1 = sigheadr1(ii).label; chan2 = sigheadr1(jj).label;
                        %        fprintf( 'iw %5d / %5d, %6s %6s \n', iw, Nw, chan1, chan2 );
                    end
                end
            end
            % AdjMatMean( :, :, iw ) = mean( AdjMatBin( :, :, 1:iw ), 3 );
        end
        
        if 1==count,
            AdjMatBin0 =  AdjMatBin ;
        else,
            lenn = length(AdjMatBin);
            lenn0 = length(AdjMatBin0);
            if lenn<lenn0,
                AdjMatBin0(:,:,lenn+1:end) = []; end;
            if lenn0<lenn,
                AdjMatBin(:,:,lenn0+1:end) = []; end;
            AdjMatBin0 = AdjMatBin0 +  AdjMatBin ;
        end
        count = count + 1;
        
    end
    
end
% Take mean
AdjMatBin = AdjMatBin0 ./ (count-1);
% initialize Adjacency matrix averaged over epochs (multiple windows, iw = 1,...,nw
AdjMatMean = zeros( size(AdjMatBin) );
for iw = 1:length(AdjMatBin),
    AdjMatMean( :, :, iw ) = mean( AdjMatBin( :, :, 1:iw ), 3 );
end

%%
if 1==iprintfig,
    
    %% PLOT ADJACENCY MATRICES
    fontsz=8;
    fh = figure(1000+ifile); set(fh,'color','w'); set(fh, 'position',[50 200 550 YLEN ]); clf; hold on;
    
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
    %colormapstrn = 'gray';
    colormap(1-gray);
    %        mymaxval = max(max(matrixplot));
    caxis([0, .1]);
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
        fnamefigconn =  ['fig_connmatrix_Tw' num2str(Tw) '_' titlestrn '.tiff'];
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefigconn] );
    end
    
    %% HEAD CONNECTIVITY
    fignum = 4000 + ifile;
    fh = figure(fignum ); set(fh,'color','w'); clf; hold on;set(fh, 'position',[600 200 350 YLEN]);
    
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
    gry = [.35 .35 .35];
    iverbose = 0; labelclr = gry; rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr );
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    drawnow
    
    if iprintfig,
        fnamefig =  ['fig_connhead_Tw' num2str(Tw) '_'  titlestrn '.tiff'];
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
    fh = figure(8000+ifile); set(fh,'color','w'); set(fh, 'position',[950 200 350 YLEN]); clf; hold on;
    
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    tccchistory =  iwhistory.*Tw;
    plot( tccchistory, ccchistory,'-','linewidth',3,'color','k');
    axis([0 Tw*Nw 0 max(ccchistory)]);
    xlabel('Time (s)','fontsize',fontsz); ylabel('Convergence of Adjacency Matrix','fontsize',fontsz);
    set(gca,'fontsize',fontsz); box on;
    
    
    if iprintfig,
        fnamefig =  ['fig_connconverge_Tw' num2str(Tw) '_'  titlestrn '.tiff'];
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
    fh = figure(9000+ifile); set(fh,'color','w'); set(fh, 'position',[1300 200 350 YLEN]); clf; hold on;
    
    % semilogy(xxxx,pppp,'.-');
    bar(xxxx,pppp, 'EdgeColor',gry,'FaceColor',gry); xlabel( 'Adjacency Weight','fontsize',fontsz ); ylabel( 'Frequency','fontsize',fontsz ); title( 'Histogram (zeros removed)' );
    if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
    axis([0 1 0 50])
    
    if iprintfig,
        fnamefig =  ['fig_connhistogram_Tw' num2str(Tw) '_' titlestrn '.tiff'];
        set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
    end
    %%
end

%% ALL FIGS IN ONE
fignum = 10+ifile;
fh = figure(fignum); set(fh,'color','w'); set(fh, 'position',[50 200 1600 YLEN ]); clf; hold on;
% PLOT ADJACENCY MATRICES
subplot('position',[0.03 0.12 .25 .74 ])

fontsz=8;
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
%colormapstrn = 'gray';
colormap(1-gray);
%        mymaxval = max(max(matrixplot));
caxis([0, .1]);
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
subplot('position',[0.34 0.12 .25 .74 ]), hold on
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
gry = [.35 .35 .35];
iverbose = 0; labelclr = gry; rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr );
if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b'); end
drawnow

% CONVERGENCE OF CONNECTIVITY
subplot('position',[0.7 0.12 .1 .74 ])

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


if 1
    % remove zero elements
    k0 = find( 0==flatvector );
    flatvector(k0) = [];
    % histogram
    [pppp,xxxx] = hist( flatvector );
    % semilogy(xxxx,pppp,'.-');
    bar(xxxx,pppp, 'EdgeColor',gry,'FaceColor',gry); xlabel( 'Adjacency Weight','fontsize',fontsz ); ylabel( 'Frequency','fontsize',fontsz ); title( 'Histogram (zeros removed)' );
    axis([0 1 0 50])
else
    [flatvecsorted, indxsorted] = sort( flatvector );
    semilogy( flatvecsorted, 'k.' );
end

if iprintfig2,
    fnamefig =  ['fig_connhistogram_Tw' num2str(Tw) '_' char(titlestrn{1}) '_' char(titlestrn{2}) '.tiff'];
    set( gcf, 'PaperPositionMode', 'auto')
    print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
end


