format compact;
fprintf('\n');
figspath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
flag_readdatafromfile = 1;
%% get data
%% compute and store
if flag_readdatafromfile,
    clear;
    dtbegeeg = 10; %sec
    dtendeeg = 10;
    filterlofreqeeg = .5;
    filterhifreqeeg = 80;
    filternotcheeg = 60;
    numsampleswindoweegpower = 200;
    filterlofreqnirs = 0.01;
    filterhifreqnirs = 0.5;
    filternotchnirs = 0;
    infostring = [ '_eegfilt' num2str(filterlofreqeeg) '-' num2str(filterhifreqeeg) '_nirsfilt' num2str(filterlofreqnirs) '-' num2str(filterhifreqnirs) ];
    infostring = [ '_eegfilt' num2str(filterlofreqeeg) '-' num2str(filterhifreqeeg) '_nirsfilt' num2str(filterlofreqnirs) '-' num2str(filterhifreqnirs) '_eegpowwindow' num2str(numsampleswindoweegpower) ];
    Fsnirs0 = 5; % Hz global sample rate to map all experiments into
    
    experimentid = [703,704,803,804,1003,1004,1102,1103,1112,1113,1122,1123,1132,1133,1134];
    
    experimentid = [1102 1112 1132 1142 1152];
    
    experimentid = [  1132 1142 1152];

    datacellsarray = cell( length(experimentid),  1 );
    for ii = 1:length(experimentid),
        % load data
        matlabdatafolder = 'D:\matlab\matlabdata';
        loadstring = ['load '  matlabdatafolder '\nirseeg_preprocessed_' num2str(experimentid(ii)) infostring '.mat'];
        fprintf('[%s]\n', loadstring);
        eval( loadstring );
        datacellsarray{ii} = datacells;
        clear datacells
    end
end

%% assign analysis params
% Global sample rate: every experiment interpolated to this before analysis
Fs0 = 5;
% Maximum time lag in delayed correlation
    maxlagtime = 20; %second

% labels of frequency bands
 freqlabel = { 'Total' '\delta' '\theta' '\alpha' '\beta' '\gamma_{tot}' '\gamma_1' '\gamma_2' };
    
% display params
flag_viewtimeseries = 0;
flag_viewtimeseries2 = 0;
flag_viewchanlocations = 0;
%%
figspath = '.\figs\';
figname = 'fig_analyzeICA_';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EEG BAND POWER - NIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over files
nirsOeventavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
nirsReventavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
nirsOblockavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
nirsRblockavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col

for experimentnum = 1:length(experimentid),
    
    datacells = datacellsarray{experimentnum};
    filenum = 1;
    %% retrieve data
    infocells = datacells{filenum,1}; infocellseeg = infocells{1}; Fs = infocellseeg{7}; eegchanlabels = infocellseeg{9}; numchaneeg = infocellseeg{8};;
    xxeeg = datacells{filenum,9}; % EEG bandpower resampled
    %xxeeg = datacells{filenum,7}; % EEG resampled
    xxeegorig = datacells{filenum,3}; % EEG original
    %xxeeg = getderivative( xxeeg );
    infocellsnirs = infocells{2}; Fsnirs = infocellsnirs{7}; nirschanlabels = infocellsnirs{9}; numchannirs = infocellsnirs{8};
    xxnirsfclip = datacells{filenum,5}; xxnirsfRclip = datacells{filenum,6}; ttnirsclip = datacells{filenum,4};
    eventtimes = datacells{filenum,10};
    %%
    badchans = infocellsnirs(6);
    %% screen display progress
    fpath = char( infocellsnirs{1} ); fpath2 = char( infocellsnirs{2} );
    fprintf('%d [%s], Fsnirs %f\n',filenum,[fpath fpath2],Fsnirs);
    %% Cconvert from MCN system to international 10-20:  T3, T4, T5 and T6— instead of T7, T8, P7 and P8
    eegchanlabels = convertMCNto1020(eegchanlabels);
    nirschanlabels = convertMCNto1020(nirschanlabels);
    %% standardize data
    
    %{
    raw:
    eegchanlabels
    nirschanlabels
    xxeeg
    xxnirsfclip
    xxnirsfRclip
    ttnirsclip
    
    standardized:
    numchans
    chanstoinclude
    eegdat
    nirsHBOdat
    nirsHBRdat
    xxxlocus
    yyylocus
    %}
    
    %%  reorder to make channels match for eeg and nirs
    
    chanstoinclude = { ...
        'FP1'...
        'FP2'...
        'F7'...
        'F8'...
        'F3'...
        'F4'...
        'FZ'...
        'T3'...
        'T4'...
        'C3'...
        'C4'...
        'CZ'...
        'T5'...
        'T6'...
        'O1'...
        'O2'...
        }; 
    chanstoinclude = { ...
        'FP1'...
        'FP2'...
        'F7'...
        'F8'...
        'F3'...
        'F4'...
        'FZ'...
        'T3'...
        'T4'...
        'C3'...
        'C4'...
        'CZ'...
        'T5'...
        'T6'...
        'P3'...
        'P4'...
        'PZ'...
        'O1'...
        'O2'...
        };
    numchans = length(chanstoinclude);
    eegindxs = [];
    nirsindxs = [];
    for kk = 1:numchans,
        mychan = chanstoinclude{kk};
        indxeeg = findstrincellarray( eegchanlabels, upper(mychan) );
        indxnirs = findstrincellarray( nirschanlabels, upper(mychan) );
        eegindxs = [eegindxs indxeeg];
        nirsindxs = [nirsindxs indxnirs];
    end
    eegdat = xxeeg(:,eegindxs,:);
    nirsHBOdat = xxnirsfclip(:,nirsindxs);
    nirsHBRdat = xxnirsfRclip(:,nirsindxs);
    
    %%  interpolate signals to a global sample rate
    TT = ttnirsclip(end);
    dt0 = 1./Fs0;
    tt0 = [0:dt0:TT];
    
    [mm,nn,kk] = size(eegdat);
    eegdat0 = zeros( length(tt0), nn, kk );
    [mm,nn] = size(nirsHBOdat);
    nirsHBOdat0 = zeros( length(tt0), nn );
    [mm,nn] = size(nirsHBRdat);
    nirsHBRdat0 = zeros( length(tt0), nn );
    
    for ichan = 1:numchans,
        
        nirsbucket = nirsHBOdat(:,ichan);
        FF = griddedInterpolant( ttnirsclip, nirsbucket, 'spline','none' );
        nirsbucket = FF(tt0)';
        nirsHBOdat0(:,ichan) = nirsbucket;
        
        nirsbucket = nirsHBRdat(:,ichan);
        FF = griddedInterpolant( ttnirsclip, nirsbucket, 'spline','none' );
        nirsbucket = FF(tt0)';
        nirsHBRdat0(:,ichan) = nirsbucket;
        
        for ifreq = 1:8,
            eegbucket = eegdat(:,ichan,ifreq);
            FF = griddedInterpolant( ttnirsclip, eegbucket, 'spline','none' );
            eegbucket = FF(tt0)';
            eegdat0(:,ichan,ifreq) = eegbucket;
        end
        if flag_viewtimeseries
            mychan = chanstoinclude{ichan};
            figure(210);clf;hold on;plot(ttnirsclip,nirsHBOdat(:,ichan),'b');plot(tt0, nirsHBOdat0(:,ichan),'r'); title(mychan,'fontsize',24);
            figure(212);clf;hold on;plot(ttnirsclip,nirsHBRdat(:,ichan),'b');plot(tt0, nirsHBRdat0(:,ichan),'r');
            figure(214);ifreq = 1; clf;hold on;plot(ttnirsclip,eegdat(:,ichan,ifreq),'b');plot(tt0,eegdat0(:,ichan,ifreq),'r');
        end
        if flag_viewtimeseries2
            figure(230);clf;hold on;
            for ifreq = 1:6,
                subplot(7,1, ifreq), plot(tt0, eegdat0(:,ichan,ifreq), 'b'); set(gca,'xtick',[])
            end
            titlestr = [ chanstoinclude{ichan} '   ' num2str(experimentid) '   ' [fpath fpath2 ] ];
            subplot(7,1,7), hold on; plot(tt0, nirsHBOdat0(:,ichan),'r');plot(tt0, nirsHBRdat0(:,ichan),'b');
            suptitle(titlestr);
        end
        
        pause(.2);
    end
        
    %% chan locations
    [ clabels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
    clabels = convertMCNto1020(clabels);
    clabels = charcellarraytoupper( clabels );
    for ii = 1:numchans,
        channam = chanstoinclude{ii};
        [ chanloctmp, iii , ifound ] =  getchanlocation( channam, clabels, xchan );
        if ifound,
            chanloc(ii,:) = chanloctmp;
        else,
            error( [channam ' not found '] );
        end
        dfas=65;
    end
    xxxlocus = chanloc(:,1); yyylocus = chanloc(:,2);
    if flag_viewchanlocations
        figure(280);clf;
        dotsiz = 22; clmap = colormap; titlestr = 'TEST'; zzz = .9.*rand(size(xxxlocus)); labelfontize = 0; fontclr = 'k';
        renderheaddots( chanstoinclude, labelfontize, fontclr, xxxlocus, yyylocus, zzz, titlestr, clmap, dotsiz, 0 );
    end
    %% calc and store delayed corr
    nlags = ceil( Fs0.*maxlagtime ); % max number of lags to compute (samples)
    maxcorrmatrixO = zeros(numchans,8);
    timemaxcorrmatrixO = zeros(numchans,8);
    maxcorrmatrixR = zeros(numchans,8);
    timemaxcorrmatrixR = zeros(numchans,8);
    
    for ichan = 1:numchans,
        nirsO = nirsHBOdat(:,ichan); nirsR = nirsHBRdat(:,ichan);
        for ifreq = 1:8,
            eegseries = eegdat(:,ichan,ifreq);
            
            [ xcorrelO, Lags ] = xcorr(  nirsO, eegseries,  nlags, 'coeff' );
            [dummy,maxindx] =  max( abs( xcorrelO ) );
            maxcorrmatrixO(ichan,ifreq) = xcorrelO(maxindx);
            timemaxcorrmatrixO(ichan,ifreq) = Lags(maxindx)./Fs0;
            
            [ xcorrelR, Lags ] = xcorr(  nirsR, eegseries,  nlags, 'coeff' );
            [dummy,maxindx] =  max( abs( xcorrelR ) );
            maxcorrmatrixR(ichan,ifreq) = xcorrelR(maxindx);
            timemaxcorrmatrixR(ichan,ifreq) = Lags(maxindx)./Fs0;
            
            if 0
                fh = figure(1020); set(fh,'color','w');clf; hold on;  grid on;
                title([ chanstoinclude{ichan} ' - ' freqlabel{ifreq} '    chan:' num2str(ichan) ' - freq:' num2str(ifreq) ],'fontsize',16);
                plot(Lags./Fs0,xcorrelO,'r'); plot( timemaxcorrmatrixO(ichan,ifreq), maxcorrmatrixO(ichan,ifreq) ,'ro' );
                plot(Lags./Fs0,xcorrelR,'b'); plot( timemaxcorrmatrixR(ichan,ifreq), maxcorrmatrixR(ichan,ifreq) ,'bo' );
                axis([-maxlagtime maxlagtime -.6 .6]);
                
                % draw channel location on head inset
                axes('position',[.7 .7 .15 .15]);
                zzz = zeros(size(xxxlocus)); zzz(ichan) = 1;
                zmax = .2; zmin = -.2; labelfontize=0; dotsiz = 10; clmap = colormap; titlestr = '';fontclr = 'k'; titlefontsz = 22; flag_axis = 0;
                renderheaddots( chanstoinclude, labelfontize, fontclr, xxxlocus, yyylocus, zzz, titlestr, titlefontsz, clmap, dotsiz, flag_axis, zmax, zmin);
                
                pause(.1);
            end
            
            corrstoreO{ichan,ifreq,experimentnum} = xcorrelO;
            corrstoreR{ichan,ifreq,experimentnum} = xcorrelR;
            
            fdas=65;
        end
        
        timelags0 = Lags./Fs0;
        
    end
    %% view maxcorr matrices
    mincorrlim = -.3;  maxcorrlim = .3;
    
    fh = figure(320); set(fh,'color','w'); clf;
    imagesc( maxcorrmatrixO );
    box off; axis off;    set(gca,'xtick',[],'ytick',[]);
    for ichan=1:numchans, text( 0, ichan, chanstoinclude{ichan} ); end
    for ifreq=1:8, text( ifreq, numchans+1, freqlabel{ifreq} ); end
    colorbar
    caxis([mincorrlim maxcorrlim]);
    
    fh = figure(330); set(fh,'color','w'); clf;
    imagesc( maxcorrmatrixR );
    box off; axis off;    set(gca,'xtick',[],'ytick',[]);
    for ichan=1:numchans, text( 0, ichan, chanstoinclude{ichan} ); end
    for ifreq=1:8, text( ifreq, numchans+1, freqlabel{ifreq} ); end
    colorbar
    caxis([mincorrlim maxcorrlim]);
    
    %% view maxcorr on the heads
    fignum = experimentid(experimentnum);
    
    zmax = .1626; zmin = -.1626;
    ysubtop = .55; ysubbot = .1; dxx = .075; dyy = .27; dotsiz = 14; clmap = colormap; titlestr = '';
    fontclr = 'k'; titlefontsz = 22; xbase = -.07; flag_axis = 0;
    
    fh = figure(fignum); set(fh,'color','w'); set(fh,'position',[100 400 1800 450]);    clf;
    for ifreq = 1:8,
        xsub = xbase + ifreq.*.115;
        axes('Position',[xsub ysubtop dxx dyy]);
        titlestr = {freqlabel{ifreq}; ' '};
        zzz = maxcorrmatrixO(:,ifreq); labelfontize=0;
        renderheaddots( chanstoinclude, labelfontize, fontclr, xxxlocus, yyylocus, zzz, titlestr, titlefontsz, clmap, dotsiz, flag_axis, zmax, zmin);
    end
    
    for ifreq = 1:8,
        xsub = xbase + ifreq.*.115; ysub = .1;
        axes('Position',[xsub ysubbot dxx dyy]);            titlestr = '';%{freqlabel{ifreq}; ' '};
        zzz = maxcorrmatrixR(:,ifreq); if ifreq==1, labelfontize=0; else, labelfontize=0; end
        renderheaddots( chanstoinclude, labelfontize, fontclr, xxxlocus, yyylocus, zzz, titlestr, titlefontsz, clmap, dotsiz, flag_axis, zmax, zmin);
    end
    % colorbar
    hb = colorbar('location','eastoutside');
    V = get(hb,'OuterPosition'); set(hb,'OuterPosition',[V(1).*1.045 V(2) V(3)./.82 V(4)./1]); set(hb,'Box','off'); caxis([zmin zmax]);
    
end

%% view full lagged corr pattern on the head
xleft = min(xxxlocus);
xright = max(xxxlocus);
ytop = max(yyylocus);
ybot = min(yyylocus);
xoffset = -.04;
width = 1.4.*(xright - xleft);
height = 1.6.*(ytop - ybot);
mydx = .12; mydy = .09;
ymaxlim = .35; yminlim = -.25;
linewdt = 1;
myxtics = [timelags0(1) 0 timelags0(end)];
myytics = [yminlim 0 ymaxlim];
fontsz = 13;

fignum = 2002;
fh = figure(fignum); set(fh,'color','w'); set(fh,'position',[100 100 800 900]);    clf;
for ichan = 1:numchans,
    chanlabl = chanstoinclude{ichan};
    myx = (xxxlocus(ichan) + width./2) ./width + xoffset; myy = (yyylocus(ichan) + height./2) ./height;
    axes( 'position', [myx myy mydx mydy] ); hold on;
    for ifreq = [6],
        
        for experimentnum = 4,%1:length(experimentid),
            curveO = corrstoreO{ichan,ifreq,experimentnum}; curveR = corrstoreR{ichan,ifreq,experimentnum};
            plot(timelags0,curveR,'b','linewidth',linewdt); plot(timelags0,curveO,'r','linewidth',linewdt);
            axis([-maxlagtime maxlagtime yminlim ymaxlim]);
            if ~isempty(findstr('O1',chanlabl)), set(gca,'XTick',myxtics); set(gca,'YTick',myytics); xlabel('Time Lag (s)'); ylabel('Correlation');
            else, set(gca,'XTick',[]); set(gca,'YTick',[]); box on;  end
            set(gca,'fontsize',fontsz); box on
            title(chanlabl);
        end
        
        pause(.5)
    end
    %    titlestr = [  num2str(experimentid) '   ' [fpath fpath2 ] ]; suptitle(titlestr);
    
end

return;
return;
return;


%{
%% delayed corr HBO-HBR
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6' 'P3' 'P4' 'PZ' }; % full 19

numchans = length(nirschanlabels);
ifreq = 6;

for ichan = 1:numchans,
    
    eegpowerseries = xxeeg(:,ichan,ifreq);
    HBOseries = xxnirsfclip(:,ichan);
    
end

%% chan-averaged EEG power for each freq band
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6' 'P3' 'P4' 'PZ' }; % full 19
flag_userelativeeegpower = true;
numfreqbands = 8;
eegpower = calc_eegpower_avgoverchans( chanstoinclude, flag_userelativeeegpower, ...
    numchaneeg, eegchanlabels, numfreqbands, xxeeg );
% view
figure(1333);clf;
for ifreqband = 1:numfreqbands,
    eegdata = eegpower{ifreqband};
    subplot(numfreqbands,1,ifreqband)
    plot(ttnirsclip, eegdata);
end
title('EEG Power time series');
%% ICA or PCA on NIRS
flag_calcICA = 0;
% select subset of NIRS chans for ICA
[nrows ncols] = size(xxnirsfclip);
indxset = [];
for ii = 1:ncols,
    channame = upper(nirschanlabels{ii}); indx = findstrincellarray( chanstoinclude, upper(channame) );
    if indx, indxset = [indxset ii]; end
end
nirschannamesubset = nirschanlabels(indxset);
numnirschansubset = length( nirschannamesubset);
xxnirsfclipsubset = xxnirsfclip(:,indxset);
xxnirsfRclipsubset = xxnirsfRclip(:,indxset);

if flag_calcICA,
    % calc IC.  x = As and s = Wx. x:signals, s: sources.
    tol = 1e-3;     % termination threshold parameter
    max_it = 1e3;   % maximum number of iterations per independent component
    % each row of W is the set of coeffs of one source
    [SS, AA, iter, WW] = robustica(xxnirsfclipsubset', [], tol, max_it, 1, 'o', 0, [], 1);   % deflationary orthogonalization
    [SSR, AAR, iterR, WWR] = robustica(xxnirsfRclipsubset', [], tol, max_it, 1, 'o', 0, [], 1);   % deflationary orthogonalization
else,
    % calc PC
    % each col of coeff contains one PC
    % each row of score is an observation in PC space
    % latent: eigenvalues
    [coeff,score,lambda,tsq] = princomp( xxnirsfclipsubset );
    WW = coeff'; SS = score';
    [coeff,score,lambda,tsq] = princomp( xxnirsfRclipsubset );
    WWR = coeff'; SSR = score';
end

%% get channel locations
[ clabels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
clabels = convertMCNto1020(clabels);
clabels = charcellarraytoupper( clabels );
for ii = 1:numnirschansubset,
    channam = nirschannamesubset{ii};
    [ chanloctmp, iii , ifound ] =  getchanlocation( channam, clabels, xchan );
    if ifound,
        chanloc(ii,:) = chanloctmp;
    else,
        error( [channam ' not found '] );
    end
    dfas=65;
end
xxx = chanloc(:,1);
yyy = chanloc(:,2);

%% view IC time series and headmaps
if 1
    for indxcomp = 1:numnirschansubset,
        
        sss = SS(indxcomp,:);
        zzz = WW(indxcomp,:)';
        
        figure(650);clf; plot( ttnirsclip, sss,'r' );
        titlestr = num2str(std(sss));
        title( titlestr, 'fontsize',16 );
        
        %% lagged corr
        fh = figure(1020); set(fh,'color','w');clf; hold on;  grid on;
        maxlagtime = 60; %second
        nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
        ifreq = 6; title([ 'comp ' num2str(indxcomp) ', freq ' num2str(ifreq) ]);
        [ xcorrel, Lags ] = xcorr(  sss, eegpower{ifreq},  nlags, 'coeff' );
        plot(Lags./Fsnirs,xcorrel,'r');
        axis([-maxlagtime maxlagtime -.6 .6]);
        %%
        
        
        fh = figure(777); set(fh,'color','w'); clf; hold on;
        %axes('Position',[.7 .7 .2 .2]);
        titlestr =  [num2str(indxcomp)]; clmap = colormap; dotsiz = 16;
        renderheaddots( xxx, yyy, zzz, titlestr, clmap, dotsiz, 0 );
        
        pause(.2);
    end
end
%% check corr of components
if 0
    for ii = 1:numnirschansubset,
        for jj = ii:numnirschansubset,
            figure(550);clf;
            myx = SSR(ii,:)'; myy = SSR(jj,:)';
            ccc = corr( [myx myy] );
            plot( myx, myy ,'b.' );
            title( [ num2str(ii) ' ' num2str(jj) ',  corr ' num2str(ccc(1,2)) ] )
            pause(.1);
        end
    end
end

asdf=65;
%% delayed corr eegpower-nirs
%xxxx = [SS' SSR' cell2mat(eegpower')'];
flag_userelativeeegpower = 1;
numfreqbands = 8;
maxlagtime = 60; %second
nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)

for indxcomp = 1:numnirschansubset,
    
    weights = abs( WW(indxcomp,:) );
    nirsseries = SS(indxcomp,:);
    zzz = WW(indxcomp,:)';
    
    %weights = ones( size( WW(indxcomp,:) ) );
    eegpower = calc_eegpower_weightedavgoverchans( chanstoinclude, flag_userelativeeegpower, ...
        numchaneeg, eegchanlabels, numfreqbands, xxeeg, weights );
    
    
    for ifreq = 1:numfreqbands,
        eegseries = eegpower{ifreq};
        [ xcorrel, Lags ] = xcorr(  nirsseries, eegseries,  nlags, 'coef' );
        
        fh = figure(1022); set(fh,'color','w');clf; hold on;  grid on;
        title([ 'comp ' num2str(indxcomp) ', freq ' num2str(ifreq) ],'fontsize',22);
        plot(Lags./Fsnirs,xcorrel,'r');
        axis([-maxlagtime maxlagtime -.6 .6]);
        
        axes('Position',[.7 .7 .2 .2]);
        titlestr =  [num2str(indxcomp)]; clmap = colormap; dotsiz = 16;
        renderheaddots( xxx, yyy, zzz, titlestr, clmap, dotsiz, 0 );
        
        pause(.1);
    end
    
end


return
%% AHMET EEGPOWER-NIRS
nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
nlags_blockavg = ceil( Fsnirs.*maxlagtime_blockavg ); % max number of lags to compute (samples)

for ii = 1:numchannirs, % NIRS
    chan1 = upper( char( nirschanlabels(ii) ) );
    
    indxbadchan = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
    if 0~=indxbadchan, continue; end
    
    indxchanstoinclude = findstrincellarray( chanstoinclude, upper(chan1) ); % return zero if chan1 is NOT a bad channel
    if 0==indxchanstoinclude, continue; end
    
    nirsOdata = ( xxnirsfclip(:,ii) );
    nirsRdata = ( xxnirsfRclip(:,ii) );
    
    for jj = 1:numchaneeg, % EEG
        chan2 = upper( char( eegchanlabels(jj) ) );
        
        if strfind( chan1, chan2 ),  % find chan2 in chan1 %  if strcmp( chan1, chan2 ),
            fprintf('[%s][%s]\n',chan1,chan2 );
            
            eegdata = ( xxeeg(:,jj,ifreqband) );
            
            if 0
                figure(9320);clf;hold on;
                [ppp,xxx] = hist(eegdata,12);
                %   ppp(1)=[];ppp(end)=[];xxx(1)=[];xxx(end)=[];
                plot( log10(xxx), log10(ppp),'.-' );
            end
            
            % choose threshold & extract times when threshold
            % was exceeded
            topfraction= .007;
            iixx = ceil(topfraction.*length(eegdata));
            [eegdatasorted, ix] = sort(eegdata,'descend');
            thresh = eegdatasorted(iixx);
            ievents = find( thresh<eegdata );
            ievents = remove_contiguous( ievents );
            
            %%
            if 1,
                figure(2321); clf; %inddxx=[1:floor(Fsnirs.*12.*60)];
                
                % eeg spectrogram
                subplot(4,1,1),
                eegdataorig = ( xxeegorig(:,jj) );
                windowspec = ceil(Fs.*1); nooverlapspec = ceil(windowspec./2);
                deltafreq = Fs./(2.*windowspec); F = [0:deltafreq:80];
                [S,F,T,P] = spectrogram(eegdataorig, windowspec, nooverlapspec, F, Fs);
                % S: Fourier coeffs, F: freqs, T: times, P: PSD
                surf(T,F,10.*log10(P),'edgecolor','none'); view(0,90); %xlabel('Time (s)');
                ylabel('Hz');
                title([chan1 ' EEG spectrogram '  fpath ' ' fpath2 ],'interpreter','none');
                ymin = 0; ymax = 80;
                %   ymin = min(eegdata); ymax = 20.*mean(eegdata);
                set(gca,'Xtick',[0:taskduration:ttnirsclip(end)]); set(gca,'TickDir','out');
                axis([ttnirsclip(1) ttnirsclip(end) ymin ymax])
                drawnow;
                
                % eeg freq band power
                subplot(4,1,2)
                gry = [.7 .7 .7];
                inddxx = [1:length(eegdata)];
                hold on; grid off;
                % event markers
                gry2 = [.7 .9 .7]; ymin = min(eegdata); ymax = max(eegdata); eventlinewidth = 2;for jjkk=1:length(eventtimes), plot( [eventtimes(jjkk) eventtimes(jjkk)], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth); end
                %plot( ttnirsclip(inddxx),xxeeg(inddxx,jj,ifreqband));
                plot(ttnirsclip(inddxx),eegdata(inddxx),'k')
                plot([ttnirsclip(1) ttnirsclip(end)],[thresh thresh],'b');
                plot(ttnirsclip(ievents),thresh.*ones(size(ievents)),'bo');
                title([ chan1 ' EEG band power '  fpath ' ' fpath2 ],'interpreter','none');
                ymin = min(eegdata);
                %   ymax = max(eegdata);
                ymax = 1000;
                %   ymax = 8.*mean(eegdata);
                set(gca,'Xtick',[0:taskduration:ttnirsclip(end)]); set(gca,'TickDir','out');
                axis([ttnirsclip(1) ttnirsclip(end) ymin ymax])
                pause(.1);
                
                %  HbO
                subplot(4,1,3)
                gry = [.7 .7 .7];
                inddxx = [1:length(eegdata)];
                hold on; grid off;
                % eeg surge markers
                plot([ttnirsclip(1) ttnirsclip(end)],[0 0],'linestyle','-','color',gry);
                ymin = min(nirsOdata); ymax = max(nirsOdata);
                %      ymin = -2e-4;  ymax = 2e-4;
                
                for jjkk=1:length(ievents),
                    %plot( [ttnirsclip(ievents(jjkk)) ttnirsclip(ievents(jjkk))], [ymin ymax],'linestyle','-','color',gry  );
                    plot( [ttnirsclip(ievents(jjkk)) ttnirsclip(ievents(jjkk))], [0 0],'linestyle','-','color',gry,'marker','o','markersize',5,'markerfacecolor',gry );
                end
                % event markers
                gry2 = [.7 .9 .7]; eventlinewidth = 2;
                for jjkk=1:length(eventtimes), plot( [eventtimes(jjkk) eventtimes(jjkk)], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth); end
                
                plot(ttnirsclip(inddxx), nirsOdata(inddxx),'r');
                set(gca,'Xtick',[0:taskduration:ttnirsclip(end)]); set(gca,'TickDir','out');
                axis([ttnirsclip(1) ttnirsclip(end) ymin ymax]);
                title([ chan1 ' HBO ' fpath ' ' fpath2 ],'interpreter','none');
                pause(.1);
                
                
                %  HbR
                subplot(4,1,4)
                gry = [.7 .7 .7];
                inddxx = [1:length(eegdata)];
                hold on; grid off;
                % event markers
                gry2 = [.7 .9 .7]; eventlinewidth = 2;
                for jjkk=1:length(eventtimes), plot( [eventtimes(jjkk) eventtimes(jjkk)], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth); end
                % eeg surge markers
                plot([ttnirsclip(1) ttnirsclip(end)],[0 0],'linestyle','-','color',gry);
                ymin = min(nirsOdata); ymax = max(nirsOdata);
                for jjkk=1:length(ievents),
                    %plot( [ttnirsclip(ievents(jjkk)) ttnirsclip(ievents(jjkk))], [ymin ymax],'linestyle','-','color',gry  );
                    plot( [ttnirsclip(ievents(jjkk)) ttnirsclip(ievents(jjkk))], [0 0],'linestyle','-','color',gry,'marker','o','markersize',5,'markerfacecolor',gry );
                end
                plot(ttnirsclip(inddxx),nirsRdata(inddxx),'b');
                set(gca,'Xtick',[0:taskduration:ttnirsclip(end)]); set(gca,'TickDir','out');
                axis([ttnirsclip(1) ttnirsclip(end) ymin ymax]);
                title([ chan1 ' HBR ' fpath ' ' fpath2 ],'interpreter','none');
                pause(.1);
                if 0,
                    fnamefig =  ['fig_EEG-NIRS_' num2str(experimentid) '.tiff'];
                    set( gcf, 'PaperPositionMode', 'auto')
                    print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
                end
                
            end
            
            %% AVERAGED NIRS SIGNAL
            % NIRS LOCKED TO EEG-SURGE-EVENT
            bucket1 = []; bucket2 = [];
            for inn = 1:length(ievents),
                indxx = [ievents(inn)-nlags:ievents(inn)+nlags];
                if 1<=indxx(1)&indxx(end)<=length(ttnirsclip),
                    bucket1 = [bucket1 nirsOdata(indxx)];
                    bucket2 = [bucket2 nirsRdata(indxx)];
                end
            end
            nirsOeventavg = mean(bucket1,2); % mean over events for one experiment one chan
            nirsReventavg = mean(bucket2,2); % mean over events for one experiment one chan
            
            if 0
                figure(1020); clf; hold on; title([ chan1 ' ' fpath ' ' fpath2 ],'interpreter','none'); grid on;
                mytim = [1:length(nirsOeventavg)]./Fsnirs-maxlagtime;
                plot( mytim, nirsOeventavg, 'r' );
                plot( mytim, nirsReventavg, 'b' );
                axis([-maxlagtime maxlagtime -2e-4 3e-4]);
                pause(.5);
            end
            
            %% RESAMPLE TO GLOBAL SAMPLE RATE (INDEPT OF EXPERIM)
            % NIRS LOCKED TO EEG-SURGE-EVENT
            mytim = [1:length(nirsOeventavg)]./Fsnirs-maxlagtime;
            nlags0 = ceil( Fsnirs0.*maxlagtime );
            mytim0 = [1:(2.*nlags0+1)]./Fsnirs0-maxlagtime;
            
            FF = griddedInterpolant( mytim, nirsOeventavg, 'spline','none' );
            nirsOeventavg0 = FF(mytim0)';
            FF = griddedInterpolant( mytim, nirsReventavg, 'spline','none' );
            nirsReventavg0 = FF(mytim0)';
            
            if 1 , % plot eeg surge event avg
                figure(1022); clf; hold on;
                title( { [ chan1]; [fpath ' ' fpath2]; ['AVG LOCKED TO EEG SURGE']},'interpreter','none','fontsize',12); grid on;
                plot( mytim0, nirsOeventavg0, 'r' );
                plot( mytim0, nirsReventavg0, 'b' );
                ymax = max([nirsOeventavg0; nirsReventavg0]); ymin = min([nirsOeventavg0; nirsReventavg0]);
                axis([-maxlagtime maxlagtime ymin ymax]);
                pause(.1);
            end
            
            nirsOeventavgbar = [nirsOeventavgbar nirsOeventavg0]; % this accumulates each chan as a column
            nirsReventavgbar = [nirsReventavgbar nirsReventavg0];% this accumulates each chan as a column
            
            
            adsf=325;
            retval = input('Press <Enter> to continue...');
            
        end
    end
end
end

%% inter experiment average
% NIRS LOCKED TO EEG-SURGE-EVENT
nirsOeventavgbar0 = mean( nirsOeventavgbar, 2 );
nirsReventavgbar0 = mean( nirsReventavgbar, 2 );
if 1
    ymin = -5e-5; ymax = 5e-5;
    
    figure(9220); clf; hold on; %title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
    plot( mytim0, nirsOeventavgbar,  'r','linewidth',1 );
    plot( mytim0, nirsOeventavgbar0,  'r','linewidth',4 );
    axis([0 maxlagtime ymin ymax]);
    
    figure(9222); clf; hold on; title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
    plot( mytim0, nirsReventavgbar,  'b','linewidth',1 );
    plot( mytim0, nirsReventavgbar0,  'b','linewidth',4 );
    axis([0 maxlagtime ymin ymax]);
    
    figure(9224); clf; hold on; title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
    plot( mytim0, nirsOeventavgbar0,  'r','linewidth',4 );
    plot( mytim0, nirsReventavgbar0,  'b','linewidth',4 );
    axis([0 maxlagtime ymin ymax]);
    
    pause(.1);
end


return;



%%
if flag_nirs_eeg_lagged_correlation,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  EEG BAND POWER - NIRS lAGGED CORR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    for kminutes = 12:12,
        
        %% loop over files
        xcorrelbar = [];%zeros(nlags.*2+1, 1);
        xcorrelRbar = [];%zeros(nlags.*2+1, 1);
        xcorrel3bar = [];%zeros(nlags.*2+1, 1);
        autocorrelEbar = [];
        autocorrelObar = [];
        autocorrelRbar = [];
        for experimentnum = [ 1 ],
            datacells = datacellsarray{experimentnum};
            filenum = 1;
            %% retrieve data
            infocells = datacells{filenum,1}; infocellseeg = infocells{1}; Fs = infocellseeg{7}; eegchanlabels = infocellseeg{9}; numchaneeg = infocellseeg{8};;
            xxeeg = datacells{filenum,9}; % EEG bandpower resampled
            %xxeeg = datacells{filenum,7}; % EEG resampled
            %xxeeg = getderivative( xxeeg );
            
            infocellsnirs = infocells{2}; Fsnirs = infocellsnirs{7}; nirschanlabels = infocellsnirs{9}; numchannirs = infocellsnirs{8};
            xxnirsfclip = datacells{filenum,5}; xxnirsfRclip = datacells{filenum,6}; ttnirsclip = datacells{filenum,4};
            %% take  part of data
            numsamptouse = floor(kminutes.*60.*Fsnirs); % floor(numsamptot./1);
            xxeeg = xxeeg(1:numsamptouse,:,:);
            xxnirsfclip = xxnirsfclip(1:numsamptouse,:);
            xxnirsfRclip = xxnirsfRclip(1:numsamptouse,:);
            %%
            badchans = infocellsnirs(6);
            %% screen display progress
            fpath = char( infocellsnirs{1} ); fpath2 = char( infocellsnirs{2} );
            fprintf('%d [%s]\n',filenum,[fpath fpath2]);
            %% Cconvert from MCN system to international 10-20:  T3, T4, T5 and T6— instead of T7, T8, P7 and P8
            eegchanlabels = convertMCNto1020(eegchanlabels);
            nirschanlabels = convertMCNto1020(nirschanlabels);
            %% Lagged corr EEGPOWER-NIRS
            nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
            
            for ii = 1:numchannirs, % NIRS
                chan1 = upper( char( nirschanlabels(ii) ) );
                
                indxbadchan = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
                if 0~=indxbadchan, continue; end
                
                indxchanstoinclude = findstrincellarray( chanstoinclude, upper(chan1) ); % return zero if chan1 is NOT a bad channel
                if 0==indxchanstoinclude, continue; end
                
                nirsOdata = meansubtract( xxnirsfclip(:,ii) );
                nirsRdata = meansubtract( xxnirsfRclip(:,ii) );
                
                [ xcorrel3, Lags ] = xcorr( nirsRdata, nirsOdata, nlags, 'coeff' );
                xcorrel3bar = [xcorrel3bar xcorrel3];
                
                [ autocorrelO, Lags ] = xcorr( nirsOdata, nirsOdata, nlags, 'coeff' );
                autocorrelObar = [autocorrelObar autocorrelO];
                
                [ autocorrelR, Lags ] = xcorr( nirsRdata, nirsRdata, nlags, 'coeff' );
                autocorrelRbar = [autocorrelRbar autocorrelR];
                
                %figure(1020);clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel3,'k');
                flagfresh = 1;
                
                for jj = 1:numchaneeg, % EEG
                    chan2 = upper ( char( eegchanlabels(jj) ) );
                    
                    if strcmp( chan1, chan2 ),
                        fprintf('[%s][%s]\n',chan1,chan2 );
                        
                        eegdata = meansubtract( xxeeg(:,jj,ifreqband) );
                        
                        if 1
                            figure(2321); clf; inddxx=[1:floor(Fsnirs.*600)]; hold on; grid on;
                            %plot( ttnirsclip(inddxx),xxeeg(inddxx,jj,ifreqband));
                            plot( ttnirsclip(inddxx), eegdata(inddxx), 'r' );
                            title([ 'EEG ' chan2 ' ' fpath ' ' fpath2 ],'interpreter','none');
                            pause(.1);
                        end
                        
                        if 0, % scatter plots
                            figure(876);clf;
                            plot( xxeeg(82:end,jj,ifreqband), xxnirsfclip(1:end-81,ii),'.');
                            title(chan2,'fontsize',24);
                            axis([-2 11 -.0003 .0003])
                            figure(8761);hold on;
                            plot(  xxnirsfRclip(92:end,ii),  xxnirsfclip(1:end-91,ii),'.'); title(chan2,'fontsize',24);
                            fsda=324;
                        end
                        
                        if 1
                            % lagged corr
                            [ xcorrel, Lags ] = xcorr(  nirsOdata, eegdata,  nlags, 'coef' );
                            [ xcorrelR, Lags ] = xcorr(  nirsRdata, eegdata, nlags, 'coef' );
                            xcorrelbar = [xcorrelbar xcorrel];
                            xcorrelRbar = [xcorrelRbar xcorrelR];
                        else
                            % lagged mutual info
                            [xcorrel, Lags] = laggedmutualinfo(  xxnirsfclip(:,ii)', xxeeg(:,jj,ifreqband)',  nlags );
                            [xcorrelR, Lags] = laggedmutualinfo(  xxnirsfRclip(:,ii)', xxeeg(:,jj,ifreqband)',  nlags );
                            xcorrelbar = [xcorrelbar xcorrel'];
                            xcorrelRbar = [xcorrelRbar xcorrelR'];
                        end
                        
                        if flagfresh,
                            % autocorr
                            [ autocorrelE, Lags ] = xcorr( eegdata, eegdata,  nlags, 'coeff' );
                            autocorrelEbar = [autocorrelEbar autocorrelE];
                        end
                        
                        if 1
                            figure(1020); clf; hold on; title([ chan1 ' ' fpath ' ' fpath2 ],'interpreter','none'); grid on;
                            plot(Lags./Fsnirs,xcorrel,'r'); plot(Lags./Fsnirs,xcorrelR,'b');
                            axis([-maxlagtime maxlagtime -.6 .6])
                            pause(.1);
                            retval = input('Press <Enter> to continue...');
                            
                        end
                    end
                end
                flagfresh = 0;
            end
        end
        xcorrelbar0 = mean(xcorrelbar,2);
        xcorrelRbar0 = mean(xcorrelRbar,2);
        xcorrel3bar0 = mean(xcorrel3bar,2);
        autocorrelEbar0 = mean(autocorrelEbar,2);
        autocorrelObar0 = mean(autocorrelObar,2);
        autocorrelRbar0 = mean(autocorrelRbar,2);
        figure(3040);clf; hold on; plot( Lags./Fsnirs, xcorrelbar ,'r' ); plot( Lags./Fsnirs, xcorrelbar0,'linewidth',5,'color','r' )
        figure(3050);clf; hold on; plot(Lags./Fsnirs, xcorrelRbar, 'b' ); plot( Lags./Fsnirs,xcorrelRbar0,'linewidth',5,'color','b' )
        %%
        fh = figure(fignum);clf; hold on;
        set(fh,'position',[100 300 600 500]); set(fh,'color','w')
        offst = -.4; ymin = -.6; ymax = .34; xmin = -maxlagtime; xmax = maxlagtime;
        gryr = [.99 .5 .5];    gryb = [.5 .5 .99]; gry = [.5 .5 .5];
        widththin = .3;
        plot([0 0],[ymin ymax],'color','k','linewidth',2,'linestyle','--'); hold on;
        plot( Lags./Fsnirs, xcorrelbar ,'color',gryr,'linewidth',widththin,'linestyle','-' );
        plot( Lags./Fsnirs, xcorrelbar0,'linewidth',5,'color','r' )
        plot(Lags./Fsnirs, xcorrelRbar + offst ,'color',gryb,'linewidth',widththin,'linestyle','-' );
        plot( Lags./Fsnirs,xcorrelRbar0+ offst,'linewidth',5,'color','b' )
        % scale bar
        plot([35 45 ],[.2 .2],'k','linewidth',3)
        %        plot([-50 -50 ],[.15 .25],'k','linewidth',3)
        text(30,.26,'10 s','fontsize',18)
        %        text(-60,.2,'0.1','fontsize',18)
        axis([xmin xmax ymin ymax])
        %set(gca,'xtick',[]);    set(gca,'ytick',[])
        set(gca,'Visible','off')
        drawnow
        if 0
            figname = 'fig_eegnirslaggedcorr_1';
            print(  '-dtiff', '-zbuffer', '-r600',  [figspath  figname  '.tiff'] );
        end
        %%
        figure(3012);clf;
        subplot(2,1,1),hold on; grid on;
        plot(Lags./Fsnirs,xcorrelbar0,'r','linewidth',3); plot(Lags./Fsnirs,xcorrelRbar0,'b','linewidth',3);
        subplot(2,1,2),hold on; grid on;
        plot(Lags./Fsnirs,autocorrelEbar0,'k','linewidth',1);
        plot(Lags./Fsnirs,autocorrelObar0,'r','linewidth',1);
        plot(Lags./Fsnirs,autocorrelRbar0,'b','linewidth',1);
        title([num2str(kminutes) ' minutes'],'fontsize',26);
        drawnow
        pause(.8)
        
        %%
        figure(3030);clf; hold on; grid on; plot(Lags./Fsnirs,xcorrel3bar0,'k');
    end
    return
end


%%
if flag_nirs_lockedto_eegfluctuations,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  EEG BAND POWER - NIRS
    %  NIRS LOCKED TO LARGE FLUCTUATIONS IN EEG POWER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    for kminutes = 14:14,
        
        %% loop over files
        nirsOeventavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
        nirsReventavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
        for experimentnum = [ 1 ],
            datacells = datacellsarray{experimentnum};
            filenum = 1;
            %% retrieve data
            infocells = datacells{filenum,1}; infocellseeg = infocells{1}; Fs = infocellseeg{7}; eegchanlabels = infocellseeg{9}; numchaneeg = infocellseeg{8};;
            xxeeg = datacells{filenum,9}; % EEG bandpower resampled
            %xxeeg = datacells{filenum,7}; % EEG resampled
            %xxeeg = getderivative( xxeeg );
            infocellsnirs = infocells{2}; Fsnirs = infocellsnirs{7}; nirschanlabels = infocellsnirs{9}; numchannirs = infocellsnirs{8};
            xxnirsfclip = datacells{filenum,5}; xxnirsfRclip = datacells{filenum,6}; ttnirsclip = datacells{filenum,4};
            %% take  part of data (to study convergence)
            numsamptouse = floor(kminutes.*60.*Fsnirs); % floor(numsamptot./1);
            numsamptouse = length( xxeeg );
            xxeeg = xxeeg(1:numsamptouse,:,:);
            xxnirsfclip = xxnirsfclip(1:numsamptouse,:);
            xxnirsfRclip = xxnirsfRclip(1:numsamptouse,:);
            ttnirsclip = ttnirsclip(1:numsamptouse);
            %%
            badchans = infocellsnirs(6);
            %% screen display progress
            fpath = char( infocellsnirs{1} ); fpath2 = char( infocellsnirs{2} );
            fprintf('%d [%s]\n',filenum,[fpath fpath2]);
            %% Cconvert from MCN system to international 10-20:  T3, T4, T5 and T6— instead of T7, T8, P7 and P8
            eegchanlabels = convertMCNto1020(eegchanlabels);
            nirschanlabels = convertMCNto1020(nirschanlabels);
            %% EMG
            figure(23211); clf; count = 1;
            for jj = 1:numchaneeg, % EMG
                chan2 = upper ( char( eegchanlabels(jj) ) );
                if [strfind( 'LOC', chan2 ) ...
                        strfind( 'CH24', chan2 ) ...
                        strfind( 'ROC', chan2 ) ...
                        strfind( 'CH25', chan2 )] ...
                        fprintf('[%s]\n',chan2 );
                    eegdata = ( xxeeg(:,jj,ifreqband) );
                    if 1,
                        inddxx=[1:floor(Fsnirs.*12.*60)]; hold on; grid off;
                        inddxx=[1:length(eegdata)];
                        subplot(4,1,count),
                        count = count + 1;
                        %plot( ttnirsclip(inddxx),xxeeg(inddxx,jj,ifreqband));
                        plot(ttnirsclip(inddxx),eegdata(inddxx),'k')
                        title([ 'EMG ' chan2 ' ' fpath ' ' fpath2 ],'interpreter','none');
                        ymin = min(eegdata); ymax = max(eegdata);
                        ymin = -1;  ymax = 60;
                        %       xmin = 0; xmax = 870;
                        axis([ttnirsclip(1) ttnirsclip(end) ymin ymax])
                        pause(.1);
                        %retval = input('Press <Enter> to continue...');
                    end
                end
            end
            
            if 0,
                fnamefig =  ['fig_EMG_' num2str(experimentid) '.tiff'];
                set( gcf, 'PaperPositionMode', 'auto')
                print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
            end
            %% EEGPOWER-NIRS
            nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
            
            for ii = 1:numchannirs, % NIRS
                chan1 = upper( char( nirschanlabels(ii) ) );
                
                indxbadchan = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
                if 0~=indxbadchan, continue; end
                
                indxchanstoinclude = findstrincellarray( chanstoinclude, upper(chan1) ); % return zero if chan1 is NOT a bad channel
                if 0==indxchanstoinclude, continue; end
                
                nirsOdata = ( xxnirsfclip(:,ii) );
                nirsRdata = ( xxnirsfRclip(:,ii) );
                
                %%
                for jj = 1:numchaneeg, % EEG
                    chan2 = upper ( char( eegchanlabels(jj) ) );
                    
                    if strfind( chan1, chan2 ),
                        fprintf('[%s][%s]\n',chan1,chan2 );
                        
                        eegdata = ( xxeeg(:,jj,ifreqband) );
                        
                        if 0
                            figure(9320);clf;hold on;
                            [ppp,xxx] = hist(eegdata,12);
                            %   ppp(1)=[];ppp(end)=[];xxx(1)=[];xxx(end)=[];
                            plot( log10(xxx), log10(ppp),'.-' );
                        end
                        
                        % choose threshold & extract times when threshold
                        % was exceeded
                        topfraction= .007;
                        iixx = ceil(topfraction.*length(eegdata));
                        [eegdatasorted, ix] = sort(eegdata,'descend');
                        thresh = eegdatasorted(iixx);
                        ievents = find( thresh<eegdata );
                        ievents = remove_contiguous( ievents );
                        
                        %    clear ievents
                        %    ievents = 4470;
                        % file  1103 event 4470 < hemo response weak although
                        % gamma strong in F. Weird doublehumpback O and P response
                        % although gamma weak there.
                        % file 1113 event 4412
                        % file 804 event 4083
                        % file 704 event 2447
                        % file 1004 event 4040;
                        % file 1004 event 3460;
                        % file 1102 ievent 3327 <- has oncenter offsurround
                        
                        %%
                        if 1,
                            figure(2321); clf;
                            
                            % eeg
                            subplot(3,1,1),
                            inddxx=[1:floor(Fsnirs.*12.*60)]; hold on; grid off;
                            inddxx=[1:length(eegdata)];
                            %plot( ttnirsclip(inddxx),xxeeg(inddxx,jj,ifreqband));
                            plot(ttnirsclip(inddxx),eegdata(inddxx),'k')
                            plot([ttnirsclip(1) ttnirsclip(end)],[thresh thresh],'b');
                            plot(ttnirsclip(ievents),thresh.*ones(size(ievents)),'bo');
                            title([  chan1 ' ' fpath ' ' fpath2 ],'interpreter','none');
                            ymin = min(eegdata); ymax = max(eegdata);
                            ymin = -1;  ymax = 60;
                            axis([ttnirsclip(1) ttnirsclip(end) ymin ymax])
                            pause(.1);
                            
                            
                            % full HbO
                            subplot(3,1,2),
                            gry = [.7 .7 .7];
                            inddxx=[1:floor(Fsnirs.*12.*60)]; hold on; grid off;
                            inddxx=[1:length(eegdata)];
                            plot([ttnirsclip(1) ttnirsclip(end)],[0 0],'linestyle','-','color',gry);
                            ymin = min(nirsOdata); ymax = max(nirsOdata);
                            xmin = ttnirsclip(1); xmax = ttnirsclip(end);
                            %ymin = -2e-4;  ymax = 2e-4;
                            %xmin = 620; xmax = 850;
                            for jjkk=1:length(ievents),
                                plot( [ttnirsclip(ievents(jjkk)) ttnirsclip(ievents(jjkk))], [ymin ymax],'linestyle','-','color',gry  );
                            end
                            plot(ttnirsclip(inddxx),nirsOdata(inddxx),'r');
                            axis([ xmin xmax ymin ymax]);
                            %title([ 'HBO ' chan1 ' ' fpath ' ' fpath2 ],'interpreter','none');
                            pause(.1);
                            
                            
                            % full Hbr
                            subplot(3,1,3),
                            gry = [.7 .7 .7];
                            inddxx=[1:floor(Fsnirs.*12.*60)]; hold on; grid off;
                            inddxx=[1:length(eegdata)];
                            plot([ttnirsclip(1) ttnirsclip(end)],[0 0],'linestyle','-','color',gry);
                            ymin = min(nirsOdata); ymax = max(nirsOdata);
                            for jjkk=1:length(ievents),
                                plot( [ttnirsclip(ievents(jjkk)) ttnirsclip(ievents(jjkk))], [ymin ymax],'linestyle','-','color',gry  );
                            end
                            plot(ttnirsclip(inddxx),nirsRdata(inddxx),'b');
                            axis([ttnirsclip(1) ttnirsclip(end) ymin ymax]);
                            %title( ['HBR ' chan1 ' ' fpath ' ' fpath2 ],'interpreter','none');
                            pause(.1);
                            
                            
                            if 0,
                                fnamefig =  ['fig_EEG-NIRS_' num2str(experimentid) '.tiff'];
                                set( gcf, 'PaperPositionMode', 'auto')
                                print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
                            end
                            
                        end
                        %%
                        
                        % fill event averaged nirs signals, nirsOeventavg
                        bucket1 = []; bucket2 = [];
                        for inn = 1:length(ievents),
                            indxx = [ievents(inn)-nlags:ievents(inn)+nlags];
                            if 1<=indxx(1)&indxx(end)<=length(ttnirsclip),
                                bucket1 = [bucket1 nirsOdata(indxx)];
                                bucket2 = [bucket2 nirsRdata(indxx)];
                            end
                        end
                        nirsOeventavg = mean(bucket1,2); % mean over events for one experiment one chan
                        nirsReventavg = mean(bucket2,2); % mean over events for one experiment one chan
                        
                        if 0
                            figure(1020); clf; hold on; title([ chan1 ' ' fpath ' ' fpath2 ],'interpreter','none'); grid on;
                            mytim = [1:length(nirsOeventavg)]./Fsnirs-maxlagtime;
                            plot( mytim, nirsOeventavg, 'r' );
                            plot( mytim, nirsReventavg, 'b' );
                            axis([-maxlagtime maxlagtime -2e-4 3e-4]);
                            pause(.5);
                        end
                        
                        % resample to global time points Fsnirs0 (indept of
                        % experiment)
                        mytim = [1:length(nirsOeventavg)]./Fsnirs-maxlagtime;
                        nlags0 = ceil( Fsnirs0.*maxlagtime );
                        mytim0 = [1:(2.*nlags0+1)]./Fsnirs0-maxlagtime;
                        
                        FF = griddedInterpolant( mytim, nirsOeventavg, 'spline','none' );
                        nirsOeventavg0 = FF(mytim0)';
                        FF = griddedInterpolant( mytim, nirsReventavg, 'spline','none' );
                        nirsReventavg0 = FF(mytim0)';
                        
                        if 1
                            figure(1022); clf; hold on; title([ chan1 ' ' fpath ' ' fpath2 ],'interpreter','none'); grid on;
                            plot( mytim0, nirsOeventavg0, 'r' );
                            plot( mytim0, nirsReventavg0, 'b' );
                            axis([-maxlagtime maxlagtime -2e-4 3e-4]);
                            title(['NIRS locked to EEG ' chan1 ' ' fpath ' ' fpath2 ],'interpreter','none');
                            pause(.5);
                            % retval = input('Press <Enter> to continue');
                        end
                        nirsOeventavgbar = [nirsOeventavgbar nirsOeventavg0]; % this accumulates each chan as a column
                        nirsReventavgbar = [nirsReventavgbar nirsReventavg0];% this accumulates each chan as a column
                        adsf=325;
                        
                    end
                end
            end
        end
        
        % inter experiment average
        nirsOeventavgbar0 = mean( nirsOeventavgbar, 2 );
        nirsReventavgbar0 = mean( nirsReventavgbar, 2 );
        
        if 1
            figure(9220); clf; hold on; %title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
            plot( mytim0, nirsOeventavgbar, 'r','linewidth',1 );
            plot( mytim0, nirsOeventavgbar0, 'r','linewidth',4 );
            
            figure(9222); clf; hold on; title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
            plot( mytim0, nirsReventavgbar, 'b','linewidth',1 );
            plot( mytim0, nirsReventavgbar0, 'b','linewidth',4 );
            
            figure(19224); clf; hold on; title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
            plot( mytim0, nirsOeventavgbar0, 'r','linewidth',4 );
            plot( mytim0, nirsReventavgbar0, 'b','linewidth',4 );
            pause(.1);
        end
        
        
        
        
    end % time
end

%%
%{
    DATA STRUCTURE
Each file is one experiment.
The file is a mat file. When the mat file is loaded, it reads a cell
array datacells which is 1 by 22. THe first dimension is always 1, higher
dimensions not used. The second dimension stores data about the experiment.
The data is as follows:

%%   data
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
    datacells{count,9} = xxeeg;
    % event times (first event time 0 since data is clipped to first event)
    datacells{count,10} = eventtimes;


    %%  info data
    infocells = cell(3,1);
    infocells{1} = infocellseeg;
    infocells{2} = infocellsnirs;
    infocellseeg = cell(8,1);
    infocellsnirs = cell(8,1);

    infocellsnirs{1} = fpath;
    infocellsnirs{2} = fpath2;
    infocellsnirs{3} = fname;
    infocellseeg{4} = fpatheeg;
    infocellseeg{5} = fnameeeg;
    infocellsnirs{6} = badchans;
    infocellseeg{7} = Fs;
    infocellsnirs{8} = Fsnirs;
    
%}
%{
                if 1==ifreqs, freqlolimit = 0; freqhilimit = 80;
            elseif 2==ifreqs, freqlolimit = 0; freqhilimit = 4;
            elseif 3==ifreqs, freqlolimit = 4; freqhilimit = 8;
            elseif 4==ifreqs, freqlolimit = 8; freqhilimit = 12;
            elseif 5==ifreqs, freqlolimit = 12; freqhilimit = 30;
            elseif 6==ifreqs, freqlolimit = 30; freqhilimit = 80;
            elseif 7==ifreqs, freqlolimit = 30; freqhilimit = 50;
            elseif 8==ifreqs, freqlolimit = 50; freqhilimit = 80;
            end
%}
%}