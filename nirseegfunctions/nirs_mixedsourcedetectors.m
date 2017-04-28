clear;
format compact;
fprintf('\n');
%% params
filterlofreqnirs = 0.01;
filterhifreqnirs = 0.5;
filternotchnirs = 0;
%% data input file
fpath = 'F:\LabData\LabData\BarbourLab\';
fpath2 = 'DrBarbourLab-SUNY-resting-state-2014-05-28_004\';
fname = 'NIRS-2014-05-28_004';
fnamelayout = [fpath '32x32motor_probeInfo.mat'];
%% src-det geometric info
load(fnamelayout);
numchans = probeInfo.probes.nChannel0;
xyc = probeInfo.probes.coords_c2;
chansrcdetindx = probeInfo.probes.index_c;
xys = probeInfo.probes.coords_s2;
xyd = probeInfo.probes.coords_d2;
numsrc = probeInfo.probes.nSource0;
chansindxleft = [1:floor(numchans./2)];
chansindxright = [floor(numchans./2)+1:numchans];
if 1
    figure(102);clf; hold on;
    xmin = 1.1.*min( [xyc(:,1); xys(:,1); xyd(:,1)] );
    xmax = 1.1.*max( [xyc(:,1); xys(:,1); xyd(:,1)] );
    ymin = 1.1.*min( [xyc(:,2); xys(:,2); xyd(:,2)] );
    ymax = 1.1.*max( [xyc(:,2); xys(:,2); xyd(:,2)] );
    for ichan = chansindxleft,
        plot( xyc(ichan,1), xyc(ichan,2), 'ko' );
        plot( xyc(ichan,1), xyc(ichan,2), 'k.' );
        isrc = chansrcdetindx(ichan,1);
        idet = chansrcdetindx(ichan,2);
        plot( xys(isrc,1),xys(isrc,2),'ro' );
        plot( xyd(idet,1),xyd(idet,2),'bo' );
        plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );
        axis([xmin xmax ymin ymax]);
    end
end

%% GET / PREPROCESS NIRS DATA
% Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
filterlofreq = filterlofreqnirs;
filterhifreq = filterhifreqnirs; %
filternotch = filternotchnirs;
% GET sample rate
Fsnirs = str2num( getfieldval( [fpath fpath2 fname '_config.txt'], 'SamplingRate', '=') );
% LOAD DATA
mydat = load( [ fpath  fpath2  fname '.wl1' ]  );
xx_wl1 = nirs_compressdata( mydat, chansrcdetindx, numsrc );
mydat = load( [ fpath  fpath2  fname '.wl2' ]  );
xx_wl2 = nirs_compressdata( mydat, chansrcdetindx, numsrc );

% convert to Hb
[xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
% preprocess
xxnirsf = mybutterandiirnotchfilters( xxo, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
xxnirsfR = mybutterandiirnotchfilters( xxr, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
%% VIEW NIRS DATA
if 0
    for ii = 1:numchans,
        chan1 = num2str(ii);
        
        figure(330);clf;
        tt = [1:length(xxnirsf)]./Fsnirs;
        subplot(2,1,1),plot( tt, xxnirsf(:, ii)); title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
        subplot(2,1,2),plot( tt, xxnirsfR(:, ii) );
        
        figure(332);clf;
        tt = [1:length(xxnirsf)]./Fsnirs;
        subplot(2,1,1),plot( tt, xx_wl1(:, ii)); title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
        subplot(2,1,2),plot( tt, xx_wl2(:, ii) );
        
        figure(440);clf;
        nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fsnirs);
        title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
        
        figure(450);clf;
        nw = 3.5; pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
        title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
        
        figure(442);clf;
        nw = 3.5; pmtm( xxo(:, ii),  nw, [], Fsnirs);
        title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
        
        figure(452);clf;
        nw = 3.5; pmtm( xxr(:, ii),  nw, [], Fsnirs);
        title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
        
        figure(460);clf;
        nw = 3.5; [Pxx,F] = pmtm( xxnirsf(:, ii),  nw, [], Fsnirs);
        logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.8 );
        plot(F,logPxx); axis([0 .8 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
        
        figure(470);clf;
        nw = 3.5; [Pxx,F] = pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
        logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.8 );
        plot(F,logPxx); axis([0 .8 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
        
        figure(494);clf;
        nw = 3.5; pmtm( xx_wl1(:, ii),  nw, [], Fsnirs);
        title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
        
        figure(496);clf;
        nw = 3.5; pmtm( xx_wl2(:, ii),  nw, [], Fsnirs);
        title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
        
        pause(.1);
    end
end

%% INTER CHAN CORRELATION
ngrid = 33;

corrmat = corr( xxnirsf(:, chansindxright) ); 
[nr,nc] = size(corrmat); corrmatHbOflat = reshape( corrmat, 1, nr.*nc );
[densitycorrO, xdensitycorrO] = hist(corrmatHbOflat,ngrid);
densitycorrO = densitycorrO./sum(densitycorrO);
figure(202);clf;hold on;
plot( xdensitycorrO, densitycorrO, 'r.-' ); grid on;

corrmat = corr( xxnirsfR(:, chansindxright) ); 
[nr,nc] = size(corrmat); corrmatHbRflat = reshape( corrmat, 1, nr.*nc );
[densitycorrR, xdensitycorrR] = hist(corrmatHbRflat,ngrid); 
densitycorrR = densitycorrR./sum(densitycorrR);
figure(204);clf;hold on;
plot( xdensitycorrR, densitycorrR, 'b.-' ); grid on;

%% LAGGED CORR GROUPS
group3 = [78 79 80];
group2 = [...    
     2,...
     3,...
     4,...
     7,...
     9,...
    10,...
    11,...
    14,...
    19,...
    20,...
    21,...
    22,...
    23,...
    24,...
    25,...
    26,...
    28,...
    30,...
    31,...
    32,...
    35,...
    36,...
    37,...
    41,...
    44,...
    45,...
    47,...
    49,...
    50,...
    52,...
    53,...
    55,...
    56,...
    57,...
    59,...
    60,...
    61,...
    62,...
    63,...
    64,...
    65,...
    66,...
    67,...
    68,...
    69,...
    71,...
    72,...
    75,...
    76,...
    83,...
    84,...
    92,...
   105,...
   109,...
   110,...
   111,...
   115,...
   119,...
   126,...
   127,...
   130,...
   136,...
   137,...
   138,...
   142,...
   143,...
   144,...
   145,...
   146,...
   147];
group1 = [     5     6    15    17    18    43    54    81    82    86    87    88    89    97    99   100   108   117   118   128   129];

%% lagged corr for each group
maxlagtime = 120;
nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
xcorrelbar = [];%zeros(nlags.*2+1, 1);
figure(1030); clf; hold on;
for ii = group1,    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrelbar = [xcorrelbar xcorrel3];    
    axis([-maxlagtime maxlagtime -1 1]);
    title([ num2str(ii) ],'interpreter','none');
    plot(Lags./Fsnirs,xcorrel3,'r');grid on;
end
title('GROUP 1');

figure(1040); clf; hold on;
for ii = group2,    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrelbar = [xcorrelbar xcorrel3];    
    axis([-maxlagtime maxlagtime -1 1]);
    title([ num2str(ii) ],'interpreter','none');
    plot(Lags./Fsnirs,xcorrel3,'b');grid on;
end
title('GROUP 2');

figure(1050); clf; hold on;
for ii = group3,    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrelbar = [xcorrelbar xcorrel3];    
    axis([-maxlagtime maxlagtime -1 1]);
    title([ num2str(ii) ],'interpreter','none');
    plot(Lags./Fsnirs,xcorrel3,'m');grid on;
end
title('GROUP 3');
%% lagged corr average of each group
maxlagtime = 120;
nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
xcorrel3bar1 = [];%zeros(nlags.*2+1, 1);
for ii = group1,    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrel3bar1 = [xcorrel3bar1 xcorrel3];    
end
xcorrel3bar1 = mean(xcorrel3bar1,2);

xcorrel3bar2 = [];%zeros(nlags.*2+1, 1);
for ii = group2,    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrel3bar2 = [xcorrel3bar2 xcorrel3];    
end
xcorrel3bar2 = mean(xcorrel3bar2,2);

xcorrel3bar3 = [];%zeros(nlags.*2+1, 1);
for ii = group3,    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrel3bar3 = [xcorrel3bar3 xcorrel3];    
end
xcorrel3bar3 = mean(xcorrel3bar3,2);

figure(2030); clf; hold on;
lnwd = 3;
plot(Lags./Fsnirs,xcorrel3bar1,'k','linewidth',lnwd);grid on;
plot(Lags./Fsnirs,xcorrel3bar2,'g','linewidth',lnwd);grid on;
%plot(Lags./Fsnirs,xcorrel3bar3,'k');grid on;
axis([-maxlagtime maxlagtime -1 1]);

%% lagged corr groups topography
figure(4040);clf; hold on;
chansiz = 12;
for ichan = [chansindxleft chansindxright],
    isrc = chansrcdetindx(ichan,1);
    idet = chansrcdetindx(ichan,2);    
    if intersect( [ichan], group1 ),
        plot( xyc(ichan,1), xyc(ichan,2), 'ko','markerfacecolor','k','markersize',chansiz );        
        plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );
    elseif intersect( [ichan], group2 ),
        plot( xyc(ichan,1), xyc(ichan,2), 'go' ,'markerfacecolor','g','markersize',chansiz);       
        plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );  
   % elseif intersect( [ichan], group3 ),
    %    plot( xyc(ichan,1), xyc(ichan,2), 'mo','markerfacecolor','m' ,'markersize',chansiz);       
     %   plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );  
    end    
    plot( xys(isrc,1),xys(isrc,2),'r.' );
    plot( xyd(idet,1),xyd(idet,2),'b.' );    
    axis([xmin xmax ymin ymax]);
end
%% PCA of Hb time series
[evec, aaa, eval, TSQUARED] = pca(xxnirsf);
evalnorm = eval./sum(eval);
figure(6010); clf;
plot(cumsum(evalnorm(1:15)),'.-'); axis([1 15 0 1]);

for imode = 1:6,
    maxval = 1.1.*max(evec(:,imode)); minval = 1.1.*min(evec(:,imode));
    fignm = ['fig_mixedsourcedet_HbR_PC' num2str(imode)];
    fh = figure(6040);clf; set(fh,'color','w')
    subplot(1,2,1), hold on;
    map = colormap; mapsz = length(map); chansiz = 9; ftsz = 18;
    for ichan = [chansindxleft chansindxright],
        val = evec(ichan, imode);
        fi = (val-minval)./(maxval-minval); mycolr = map( 1+floor( fi.*mapsz ), : );
        plot( xyc(ichan,1), xyc(ichan,2), 'o','markerfacecolor',mycolr,'markeredgecolor',mycolr,'markersize',chansiz );
    end
    set(gca, 'visible', 'off')
    axis([-.8 .8 -.8 .4]); 
    text(-.8,-.8,['HbR PC ' num2str(imode) ', \lambda=' num2str(evalnorm(imode))],'fontsize',ftsz);
    
    subplot(1,2,2)
    nw = 3.5; [Pxx,F] = pmtm( aaa(:, imode),  nw, [], Fsnirs);
    logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.5 );
    plot(F,logPxx); axis([0 .5 min(logPxx(indxs)) max(logPxx(indxs))]);grid on; xlabel('Frequency (Hz)','fontsize',ftsz); 
    ylabel('Power/Freq. (dB/Hz)','fontsize',ftsz);
    
    drawnow
    if 0
        figspath = '.\figs\';                set( gcf, 'PaperPositionMode', 'auto')
        print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fignm '.tiff'] );
    end
    pause(1)
end
%% cluster lagged corr
% calc set of lagged corrs
maxlagtime = 60; numclust = 4;
nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
xcorrelbar = [];%zeros(nlags.*2+1, 1);
for ii = [chansindxleft chansindxright],    
    [ xcorrel, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrelbar = [xcorrelbar xcorrel];    
end
% cluster them
IDX = kmeans(xcorrelbar',numclust);
% cluster corr plots

fh = figure(7040); set(fh,'color','w'); clf; hold on; lnwd = 1; fnsz = 18;
for iplot = 1:numclust,
    subplot(numclust,2,iplot.*2-1), hold on;
    iclust = iplot; ind = find(iclust==IDX);
    for ii = 1:length(ind), jj = ind(ii); plot(Lags./Fsnirs,xcorrelbar(:,jj),'b','linewidth',lnwd);grid on;end
    set(gca,'XTick',[-60 0 60]);    set(gca,'YTick',[]); box on;

    if iplot==numclust,xlabel('Lag (s)','fontsize',fnsz);end
    axis([-maxlagtime maxlagtime -1 1]);
    
    gry = [.8 .8 .8]; chanclr = 'b';
    subplot(numclust,2,iplot.*2), hold on;
    chansiz = 8;
    for ichan = [chansindxleft chansindxright],
        plot( xyc(ichan,1), xyc(ichan,2), 'o','color',gry,'markerfacecolor',gry,'markersize',5 );
    end
    for ichan = [chansindxleft chansindxright],
        if IDX(ichan)==(iplot),
            plot( xyc(ichan,1), xyc(ichan,2), 'o','color',chanclr,'markerfacecolor',chanclr,'markersize',chansiz );
        end
    end
    set(gca, 'visible', 'off')
    axis([-.8 .8 -.8 .4]);
end

%%
clustcolr = ['r' 'b' 'g' 'k'];
fh = figure(7050); set(fh,'color','w'); clf; hold on; lnwd = .4; fnsz = 18;
for iplot = 1:numclust,
    subplot(numclust,1,iplot), hold on;
    iclust = iplot; ind = find(iclust==IDX);
    mycolr = clustcolr(iplot);
    for ii = 1:length(ind), jj = ind(ii); plot(Lags./Fsnirs,xcorrelbar(:,jj),mycolr,'linewidth',lnwd);grid on;end
    set(gca,'XTick',[-60 -30 0 30 60]);    set(gca,'YTick',[-1 0 1]); box on;
    if iplot==numclust,xlabel('Lag (s)','fontsize',fnsz);end
    axis([-maxlagtime maxlagtime -1 1]);   
end
drawnow
if 0
    fignm = 'fig_mixedsrcdet_clusteredcorrs';
    figspath = '.\figs\';                set( gcf, 'PaperPositionMode', 'auto')
    print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fignm '.tiff'] );
end
%%
fh = figure(7060); set(fh,'color','w'); clf; hold on; lnwd = 1; fnsz = 18;
gry = [.8 .8 .8];
chansiz1 = 5; chansiz2 = 14;
for ichan = [chansindxleft chansindxright],
    plot( xyc(ichan,1), xyc(ichan,2), 'o','color',gry,'markerfacecolor',gry,'markersize',chansiz1 );
    
    isrc = chansrcdetindx(ichan,1);
    idet = chansrcdetindx(ichan,2);
    plot( xys(isrc,1),xys(isrc,2),'ro' );
    plot( xyd(idet,1),xyd(idet,2),'bo' );
    plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );
end
for ichan = [chansindxleft chansindxright],
    for iclust = 1:numclust,
        mycolr = clustcolr(iclust);
        if iclust==IDX(ichan),
            plot( xyc(ichan,1), xyc(ichan,2), 'o','color',mycolr,'markerfacecolor',mycolr,'markersize',chansiz2 );
        end
    end
end
set(gca, 'visible', 'off')
axis([-.8 .8 -.8 .4]);

drawnow
if 0
    fignm = 'fig_mixedsrcdet_clusteredcorrs_topo';
    figspath = '.\figs\';                set( gcf, 'PaperPositionMode', 'auto')
    print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fignm '.tiff'] );
end
return;
%% topo plot of clusters
figure(7040);clf; hold on;
chansiz = 12;
for ichan = [chansindxleft chansindxright],
    isrc = chansrcdetindx(ichan,1);
    idet = chansrcdetindx(ichan,2);    
    if 1==IDX(ichan),
        plot( xyc(ichan,1), xyc(ichan,2), 'ko','markerfacecolor','k','markersize',chansiz );        
        plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );
    elseif 2==IDX(ichan),
        plot( xyc(ichan,1), xyc(ichan,2), 'go' ,'markerfacecolor','g','markersize',chansiz);       
        plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );  
   % elseif intersect( [ichan], group3 ),
    %    plot( xyc(ichan,1), xyc(ichan,2), 'mo','markerfacecolor','m' ,'markersize',chansiz);       
     %   plot( [xys(isrc,1) xyd(idet,1)] , [xys(isrc,2) xyd(idet,2)],'k:' );  
    end    
    plot( xys(isrc,1),xys(isrc,2),'r.' );
    plot( xyd(idet,1),xyd(idet,2),'b.' );    
    axis([xmin xmax ymin ymax]);
end
%% PCA of lagged corr
% calc set of lagged corrs
maxlagtime = 60;
nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
xcorrelbar = [];%zeros(nlags.*2+1, 1);
for ii = [chansindxleft chansindxright],    
    [ xcorrel, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrelbar = [xcorrelbar xcorrel];    
end
[evec, aa, eval, TSQUARED, evalcum] = pca(xcorrelbar');
evalnorm = eval./sum(eval);
figure(501); clf;
plot(cumsum(evalnorm(1:15)),'.-'); axis([1 15 0 1]);
figure(502);clf; hold on;
plot(Lags./Fsnirs,evec(:,1),'r');grid on;
plot(Lags./Fsnirs,evec(:,2),'b');grid on;
plot(Lags./Fsnirs,evec(:,3),'k');grid on;
plot(Lags./Fsnirs,evec(:,4),'g');grid on;

%%
return
%%
for ii = [chansindxleft chansindxright],
    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrel3bar = [xcorrel3bar xcorrel3];
    
    axis([-maxlagtime maxlagtime -1 1]);
    title([ num2str(ii) ],'interpreter','none');
    
    if intersect( [ii], group1 ),
        plot(Lags./Fsnirs,xcorrel3,'r');grid on;
    elseif intersect( [ii], group2 ),
        plot(Lags./Fsnirs,xcorrel3,'b');grid on;
    else
        plot(Lags./Fsnirs,xcorrel3,'g');grid on;
    end
    
end
%%
return
%% COMPUTE lagged corr GROUPS
if 0
subindx = [];subindx2 = [];subindx3 = [];

for ii = [chansindxleft chansindxright],
    
    [ xcorrel3, Lags ] = xcorr(   xxnirsfR(:,ii), xxnirsf(:,ii), nlags, 'coeff' );
    xcorrel3bar = [xcorrel3bar xcorrel3];
    plot(Lags./Fsnirs,xcorrel3,'k');grid on;
    
    axis([-maxlagtime maxlagtime -1 1]);
    title([ num2str(ii) ],'interpreter','none');
    
    indx0 = find( abs(Lags./Fsnirs)<2 );
    ind = find( (xcorrel3(indx0))>.4 );
    
    if ~isempty(ind),
        subindx = [subindx; ii];
        hold on;
        plot(Lags(indx0)./Fsnirs,xcorrel3(indx0),'ro-'); grid on;
        pause(.5)
        hold off;
        continue
    end
    
    indx0 = find( abs(Lags./Fsnirs)<2 );
    ind = find( (xcorrel3(indx0))<-.7 );
    
    if ~isempty(ind),
        subindx3 = [subindx3; ii];
        hold on;
        plot(Lags(indx0)./Fsnirs,xcorrel3(indx0),'mo-'); grid on;
        pause(.5)
        hold off;
        continue
    end
    
    indx0 = find( 4<(Lags./Fsnirs)&(Lags./Fsnirs)<18 );
    ind = find( (xcorrel3(indx0))<-.45 );
    if ~isempty(ind),
        subindx2 = [subindx2; ii];
        hold on;
        plot(Lags(indx0)./Fsnirs,xcorrel3(indx0),'bo-'); grid on;
        pause(.5)
        hold off;
    end
    
    
    
    pause(.25)
    
end
end
%%
return;

%%
for icycle = experimentid,
    
    if 703==icycle,
    elseif 1004==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            % fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\Lu0609\NIRS\';
        fpath2 = 'Resting State\';
        fname = 'NIRS-2014-06-09_004';
        fpatheeg = 'F:\LabData\LabData\Lu0609\EEG\';
        fnameeeg = 'lurestingstate_20140609_175241';
        badchans = {'CZ' };
        nirschanreplace = [];
        %{
        %}
    end
    
    fprintf('[%s]\n',[fpath fpath2 fname]  );
    fprintf('[%s]\n', [fpatheeg fnameeeg] );
    %% get EEG
    dtbeg = dtbegeeg;
    dtend = dtendeeg;
    %    dtbeg = 10; % sec
    %    dtend = 10; % sec
    % EEG  Bandpas filter and notch (normal bandpass range: 0.5-70 Hz)
    filterlofreq = filterlofreqeeg;
    filterhifreq = filterhifreqeeg; %
    filternotch = filternotcheeg;
    %    filterlofreq = 0; % .5; % Hz
    %    filterhifreq = 70; % 70; % Hz
    %    filternotch = 60; % 0: no notch;
    file1 = [fpatheeg fnameeeg '.edf'];
    [eegchanlabels, eegchanlocs, xxeegf, evteeg, Fs] = getpreprocess_eeg( file1, dtbeg, dtend, filterlofreq, filterhifreq, filternotch );
    numchaneeg = length(eegchanlabels);
    %% get NIRS
    %% GET / PREPROCESS NIRS DATA
    % Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
    filterlofreq = filterlofreqnirs;
    filterhifreq = filterhifreqnirs; %
    filternotch = filternotchnirs;
    %    filterlofreq = .01; % Hz
    %    filterhifreq = .08; % Hz
    %    filternotch = 0; % 0: no notch;[
    % GET sample rate
    Fsnirs = str2num( getfieldval( [fpath fpath2 fname '_config.txt'], 'SamplingRate', '=') );
    %  EVENTS filename
    filenirsevt = [ fpath  fpath2  fname  '.evt' ];
    
    % LOAD DATA
    mydat = load( [ fpath  fpath2  fname '.wl1' ]  );
    xx = fnnirs_reduce_date( mydat );
    xx_wl1 = xx;
    mydat = load( [ fpath  fpath2  fname '.wl2' ]  );
    xx = fnnirs_reduce_date( mydat );
    xx_wl2 = xx;
    
    % convert to Hb
    [xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
    % preprocess
    xxnirsf = mybutterandiirnotchfilters( xxo, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
    xxnirsfR = mybutterandiirnotchfilters( xxr, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
    %% VIEW NIRS DATA
    if 0
        for ii = 1:numchan,
            chan1 = upper( nirschanlabels(ii) );
            
            figure(330);clf;
            tt = [1:length(xxnirsf)]./Fsnirs;
            subplot(2,1,1),plot( tt, xxnirsf(:, ii)); title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
            subplot(2,1,2),plot( tt, xxnirsfR(:, ii) );
            
            figure(332);clf;
            tt = [1:length(xxnirsf)]./Fsnirs;
            subplot(2,1,1),plot( tt, xx_wl1(:, ii)); title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
            subplot(2,1,2),plot( tt, xx_wl2(:, ii) );
            
            figure(440);clf;
            nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            
            figure(450);clf;
            nw = 3.5; pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            
            figure(442);clf; 
            nw = 3.5; pmtm( xxo(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            
            figure(452);clf;
            nw = 3.5; pmtm( xxr(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            
            figure(460);clf;
            nw = 3.5; [Pxx,F] = pmtm( xxnirsf(:, ii),  nw, [], Fsnirs);
            logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.8 );
            plot(F,logPxx); axis([0 .8 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
                        
            figure(470);clf;
            nw = 3.5; [Pxx,F] = pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
            logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.8 );
            plot(F,logPxx); axis([0 .8 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
               
            figure(494);clf; 
            nw = 3.5; pmtm( xx_wl1(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            
            figure(496);clf; 
            nw = 3.5; pmtm( xx_wl2(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' get_chanlabel_from_srcdetpair( ii, ii )  ] , 'fontsize', 22);
            
            pause(.1);
        end
    end
    %% SYNC
    % EEG event sample numbers
    fprintf('Clipping data to sync\n');
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
        ttnirsclip = ttnirsclip(1:nnn);
        xxnirsfclip = xxnirsfclip(1:nnn, :);
        xxnirsfRclip = xxnirsfRclip(1:nnn, :);
    else
        % the other way around
        nnn = min( find( ttnirsclip(end)<tteegclip ) ) - 1;
        tteegclip = tteegclip(1:nnn);
        xxeegfclip = xxeegfclip(1:nnn, :);
    end
    %% remove NIRS freq from EEG
    if 0
        for ii = 1:numchaneeg,
            sss = xxeegfclip(:,ii);
            SSS = fft(sss);
            %nw = 3.5;
            %figure(12345); clf; pmtm( sss, nw, [], Fs);
            
            for klm = 1:14,
                nnn = ceil( klm.*length(sss).*Fsnirs./Fs );
                for nnn1=[nnn-22:nnn+22], SSS(nnn1) = 0; SSS(end-nnn1+1) = 0; end;
            end
            
            sssr = real( ifft(SSS) );
            %figure(22345); clf; pmtm( sssr, nw, [], Fs);
            
            xxeegfclip(:,ii) = sssr;
        end
    end
    
    %% INTERPOLATE EEG TO NIRS DATA POINTS
    xxeegfclipds = zeros( length(ttnirsclip), numchan );
    for ii = 1:numchaneeg,
        FF = griddedInterpolant( tteegclip, xxeegfclip(:,ii), 'spline','none' );
        xxeegfclipds(:,ii) = FF(ttnirsclip);
    end
    %% view eeg spectrum
     if 0
        for ii = 1:numchaneeg,
            chan1 = upper( eegchanlabels(ii) );
            nw=3.5;
            figure(2330);clf;
%           plot( tteegclip, xxeegfclip(:,ii) ); 
            pmtm( xxeegfclip(:,ii),nw,[],Fs);title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
            pause(1);
            
            sss = xxeegfclip(:,ii);
            SSS = fft(sss);
            figure(1234); clf; plot( log(abs(SSS)) );
            
            for klm = 1:14,
                nnn = ceil( klm.*length(sss).*Fsnirs./Fs );
                for nnn1=[nnn-12:nnn+12], SSS(nnn1) = 0; SSS(end-nnn1+1) = 0; end;
            end
            figure(12341); clf; plot( log(abs(SSS)) );
            
            sssr = real( ifft(SSS) );
                        
            figure(2345); clf; plot( tteegclip(1:2000),sssr(1:2000), tteegclip(1:2000),sss(1:2000),'r' )

            SSSr = fft(sssr);            
            figure(12345); clf; plot( log(abs(SSSr)) );
            figure(22345); clf; pmtm( sssr, nw, [], Fs);
            pause(.1);
        end
     end
    %% EEG POWER
    % calc power time course
    nweegpower = numsampleswindoweegpower;
%    nweegpower = 50; % (must be even) window size in num samples
    num = nweegpower-1; denom = nweegpower; noverlap = floor( nweegpower.* (num./denom) );
    
    numfreqbands = 6;
    xxeegfclipbandpower = zeros( length(tteegclip), numchaneeg, numfreqbands );
    xxeegfclipbandpowerds = zeros( length(ttnirsclip), numchaneeg, numfreqbands );
    
    for kchan = 1:numchaneeg,
        fprintf('Calculating eeg power, chan %d/%d\n',kchan,numchaneeg);
        freqs = [];
        [SS,Freq,Tim,Pow] = spectrogram(xxeegfclip(:,kchan),nweegpower,noverlap,freqs,Fs);
        %Pow = 10.*log10(Pow); % converting to dB has little effect on lagged corr
        
        % Pow is shorter than the original time series because of window
        % size. Replicate first/last entry to pad the missing initial/final
        % segments
        avecinitial = Pow(:,1); avecfinal = Pow(:,end);
        Pow = [ repmat(avecinitial,1,nweegpower./2-1) Pow repmat(avecfinal,1,nweegpower./2) ];
        
        % loop over frequency bands
        for ifreqs = [1:numfreqbands],
            if 1==ifreqs, freqlolimit = 0; freqhilimit = 70;
            elseif 2==ifreqs, freqlolimit = 0; freqhilimit = 4;
            elseif 3==ifreqs, freqlolimit = 4; freqhilimit = 8;
            elseif 4==ifreqs, freqlolimit = 8; freqhilimit = 12;
            elseif 5==ifreqs, freqlolimit = 12; freqhilimit = 30;
            elseif 6==ifreqs, freqlolimit = 30; freqhilimit = 70; end
            
            % remove nirs freqs
            if 0
                for klm = 1:22,
                    f1 = klm.*Fsnirs-.4; f2 = klm.*Fsnirs+.4;
                    indxs = find( f1<=Freq & Freq<=f2 );
                    Pow(indxs) = 0;
                end
            end
            
            %EEG power  ( time x chans x freqbands )
            indxs = find( freqlolimit<=Freq & Freq<=freqhilimit );
            totalbandpower = sum( Pow(indxs,:), 1 );
            xxeegfclipbandpower(:,kchan,ifreqs) = totalbandpower;
            
            %EEG power resampled to NIRS time points ( time x chans x freqbands )
            FF = griddedInterpolant( tteegclip, xxeegfclipbandpower(:,kchan,ifreqs), 'spline','none' );
            xxeegfclipbandpowerds(:,kchan,ifreqs) = FF(ttnirsclip);
        end
    end
    
    %% assign info data
    infocellseeg = cell(22,1);
    infocellsnirs = cell(22,1);
    infocells = cell(3,1);
    
    infocellsnirs{1} = fpath;
    infocellsnirs{2} = fpath2;
    infocellsnirs{3} = fname;  
    infocellsnirs{4} = badchans;
    infocellsnirs{7} = Fsnirs;
    infocellsnirs{8} = numchan;
    infocellsnirs{9} = nirschanlabels;
    
    infocellseeg{1} = fpatheeg;
    infocellseeg{3} = fnameeeg;
    infocellseeg{7} = Fs;
    infocellseeg{8} = numchaneeg;
    infocellseeg{9} = eegchanlabels;
    
    infocells{1} = infocellseeg;
    infocells{2} = infocellsnirs;
    %% assign return data
    datacells{count,1} = infocells;
    % synced data
    datacells{count,2} = tteegclip;
    datacells{count,3} = xxeegfclip;
    datacells{count,4} = ttnirsclip;
    datacells{count,5} = xxnirsfclip;
    datacells{count,6} = xxnirsfRclip;
    % eeg resampled to nirs time points
    datacells{count,7} = xxeegfclipds;
    % eeg band power
    datacells{count,8} = xxeegfclipbandpower;
    % eeg band power resampled to nirs time points
    datacells{count,9} = xxeegfclipbandpowerds;
        
    count = count + 1;
    
end

