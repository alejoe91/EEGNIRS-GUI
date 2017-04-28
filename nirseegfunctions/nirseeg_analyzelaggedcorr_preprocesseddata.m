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
    experimentid = [    1142  ];
    
    datacellsarray = cell( length(experimentid),  1 );
    for ii = 1:length(experimentid),
        % load data
        matlabdatafolder = '\matlab\matlabdata';
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
maxlagtime = 45; %second

% labels of frequency bands
freqlabel = { 'Total' '\delta' '\theta' '\alpha' '\beta' '\gamma_{tot}' '\gamma_1' '\gamma_2' };

% display params
flag_viewtimeseries = 0;
flag_viewtimeseries2 = 0;
flag_viewchanlocations = 0;
flag_viewindividuallaggedcorr = 0;
flag_viewmaximumcorrmatrix = 0;
flag_viewmaxcorronhead = 0;

% channels to use
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
%%
figspath = '.\figs\';
figname = 'fig_analyzelaggedcorr_';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EEG BAND POWER - NIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over files
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
        pause(.8);
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
            
            if flag_viewindividuallaggedcorr,
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
                
                pause(.8);
            end
            
            corrstoreO{ichan,ifreq,experimentnum} = xcorrelO;
            corrstoreR{ichan,ifreq,experimentnum} = xcorrelR;
            
        end
        
        timelags0 = Lags./Fs0;
        
    end
    %% view maxcorr matrices
    if flag_viewmaximumcorrmatrix,
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
    end
    
    %% view maxcorr on the heads
    if flag_viewmaxcorronhead,
        
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
fontsz = 12;
fontsztitle = 16;

fignum = 2002;
fh = figure(fignum); set(fh,'color','w'); set(fh,'position',[100 100 800 900]);    clf;
for ichan = 1:numchans,
    chanlabl = chanstoinclude{ichan};
    myx = (xxxlocus(ichan) + width./2) ./width + xoffset; myy = (yyylocus(ichan) + height./2) ./height;
    axes( 'position', [myx myy mydx mydy] ); hold on;
    for ifreq = [6],
        
        for experimentnum = 1,%1:length(experimentid),
            curveO = corrstoreO{ichan,ifreq,experimentnum}; curveR = corrstoreR{ichan,ifreq,experimentnum};
            plot(timelags0,curveR,'b','linewidth',linewdt); plot(timelags0,curveO,'r','linewidth',linewdt);
            axis([-maxlagtime maxlagtime yminlim ymaxlim]);
            if ~isempty(findstr('O1',chanlabl)), 
                set(gca,'XTick',myxtics,'fontsize',fontsz); set(gca,'YTick',myytics,'fontsize',fontsz); 
                xlabel('Time Lag (s)'); ylabel('Correlation');
            else, set(gca,'XTick',[]); set(gca,'YTick',[]); box on;  end
            set(gca,'fontsize',fontsz); box on
            title(chanlabl,'fontsize',fontsztitle);
        end
        
        pause(.01)
    end
    %    titlestr = [  num2str(experimentid) '   ' [fpath fpath2 ] ]; suptitle(titlestr);
    
end

return;


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