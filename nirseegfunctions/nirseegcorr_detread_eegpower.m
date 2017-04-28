clear; fprintf('\n'); format compact;
%% PARAMS
fignum = 42311;
idisplaytitle = 1;
Tw = 40; % 40; % analysis window (sec) 40
maxlagtime = 30;
flag_print_fig = 0;
figspath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
%figname = 'fig_eeg_nirs_laggedcorr_alpha_freq8-12_task_301';
%figname = 'fig_eeg_nirs_laggedcorr_broadband_task_301';
%figname = 'fig_eeg_nirs_laggedcorr_beta_freq16-31_task_301';
figname = 'fig_eegpower_';
flag_useprefrontalonly = 0;

%% EEG POWER PARAMS
freqlolimit = .5; freqhilimit = 70;
nweegpower = 50; % window size in num samples
%% choose data file
%%
%%
for icycle = [1004],
    
    %% EEG  Bandpas filter and notch (normal bandpass range: 0.5-70 Hz)
    filterlofreq = 0; % .5; % Hz
    filterhifreq = 70; % 70; % Hz
    filternotch = 60; % 0: no notch;
    flagytics = 0; % disp y lables and tics
    flagxtics = 0; % disp x labels and tics
    % axis limits
    xminlim = -maxlagtime; xmaxlim = maxlagtime; yminlim = -.01; ymaxlim = .02;
    xminlim = -22; xmaxlim = 22; yminlim = -.01; ymaxlim = .018;
    xminlim = -25; xmaxlim = 25; yminlim = -.2; ymaxlim = .2;
    % fontsizes
    fontsz = 26; fontsztic = 16; linsz = 2; linsz2 = 2;
    %% DATA FILES
    if 9==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\RestingState\';
        fname = 'NIRS-2014-03-27_002';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
        fnameeeg = ['Omurtag327restingstateboth_20140327_162652'];
        badchans = {}; % example: badchans = {'FP1', 'T6' 'F4' 'FP1' 'F1' 'F8' 'P3' 'EKG1' };
        nirschanreplace = [17 20];
    elseif 10==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\RestingStateeyesopen\';
        fname = 'NIRS-2014-03-27_003';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
        fnameeeg = ['Omurtag327restingstateeyesopenboth_20140327_163413'];
        nirschanreplace = [17 20];
        badchans = {};
    elseif 11==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\Taskeyesclose\';
        fname = 'NIRS-2014-03-27_004';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
        fnameeeg = ['Omurtag327taskeyescloseboth_20140327_164238'];
        nirschanreplace = [17 20];
        badchans = {};
    elseif 12==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0327\NIRS\';
        fpath2 = 'Detector Readings\Taskeyesopen\';
        fname = 'NIRS-2014-03-27_005';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
        fnameeeg = ['Omurtag327taskeyesopenboth_20140327_165145'];
        badchans = {};
        nirschanreplace = [17 20];
    elseif 201==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0415\NIRS\';
        fpath2 = 'Detector Readings\Task\';
        fname = 'NIRS-2014-04-14_001';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0415\EEG\';
        fnameeeg = ['AhmetOmurtagtask15min_20140414_185821'];
        badchans = {};
        nirschanreplace = [17 20];
    elseif 301==icycle,
        numchan = 19;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0416\NIRS\';
        fpath2 = 'Detector Readings\task\';
        fname = 'NIRS-2014-04-16_003';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0416\EEG\';
        fnameeeg = ['AhmetOmurtag41615minutestask_20140416_191933'];
        badchans = {};
        nirschanreplace = [17 20];
    elseif 302==icycle,
        numchan = 19;
        for ii = 1:numchan,        nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );        fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0416\NIRS\';
        fpath2 = 'Detector Readings\eyesclose\';
        fname = 'NIRS-2014-04-16_005';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0416\EEG\';
        fnameeeg = ['AhmetOmurtag41615minuteseyesclosefinal_20140416_200010'];
        badchans = {};
        nirschanreplace = [17 20];
    elseif 303==icycle,
        numchan = 19;
        for ii = 1:numchan,        nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );        fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0416\NIRS\';
        fpath2 = 'Detector Readings\eyesopen\';
        fname = 'NIRS-2014-04-16_004';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0416\EEG\';
        fnameeeg = ['AhmetOmurtag41615minuteseyesopenfinal_20140416_194254'];
        badchans = {};
        nirschanreplace = [17 20];
    elseif 401==icycle,
        numchan = 19;
        for ii = 1:numchan,        nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );        fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\Thinh Nguyen0418\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-04-18_001';
        fpatheeg = 'F:\LabData\LabData\Thinh Nguyen0418\EEG\';
        fnameeeg = ['thintask10min_20140418_154457'];
        badchans = {};
        nirschanreplace = [17 20];
        %{
    Kullanilabilir kanallar : 1-2-3-4-9-16-18-19
Kullanilabilme ihtimali olanlar ( cok gurultulu) - 6-7-8-14-15-16-10
Kullanimi imkansiz: 5-11-12-13-17
        %}
    elseif 402==icycle,
        numchan = 19;
        for ii = 1:numchan,        nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );        fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\Thinh Nguyen0418\NIRS\';
        fpath2 = 'Resting\';
        fname = 'NIRS-2014-04-18_003';
        fpatheeg = 'F:\LabData\LabData\Thinh Nguyen0418\EEG\';
        fnameeeg = ['thinresting10min_20140418_155803'];
        badchans = {};
        nirschanreplace = [17 20];
        %{
    Kullanilabilir kanallar : 1-2-3-4-9-16-18-19
Kullanilabilme ihtimali olanlar ( cok gurultulu) - 6-7-8-14-15-16-10
Kullanimi imkansiz: 5-11-12-13-17
        %}
    elseif 501==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\ErgunCumbul0421\NIRS\';
        fpath2 = 'task\';
        fname = 'NIRS-2014-04-20_004';
        fpatheeg = 'F:\LabData\LabData\ErgunCumbul0421\EEG\';
        fnameeeg = ['erguntask15minfinalfinal_20140420_215443'];
        badchans = {};
        nirschanreplace = [17 20; 16 21];
        %{
SD17 yerine SD20
SD16 yerine SD21
        %}
    elseif 502==icycle,
        numchan = 19;
        for ii = 1:numchan, chanlabel = get_chanlabel_from_srcdetpair( ii, ii ); fprintf( '%d %d  %s\n', ii, ii, chanlabel );end
        fpath = 'F:\LabData\LabData\ErgunCumbul0421\NIRS\';
        fpath2 = 'resting state\';
        fname = 'NIRS-2014-04-20_005';
        fpatheeg = 'F:\LabData\LabData\ErgunCumbul0421\EEG\';
        fnameeeg = ['ergunresting15min_20140420_221123'];
        badchans = {};
        nirschanreplace = [17 20; 16 21];
        %{
SD17 yerine SD20
SD16 yerine SD21
        %}
    elseif 601==icycle,
        numchan = 19;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            %        nirschanlabels{ii} =         get_chanlabel_from_srcdec(ii,ii);
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0509\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-05-09_003';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0509\EEG\';
        fnameeeg = 'aomurtagtask_20140509_175335';
        badchans = {'F8' 'T8' 'P8'};
        nirschanreplace = [];
    elseif 602==icycle,
        numchan = 19;
        for ii = 1:numchan,        nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );        fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        
        fpath = 'F:\LabData\LabData\AhmetOmurtag0509\NIRS\';
        fpath2 = 'Resting\';
        fname = 'NIRS-2014-05-09_004';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0509\EEG\';
        fnameeeg = 'aomurtagrestingp_20140509_181044';
        badchans = {'F8' 'T8' 'P8'};
        nirschanreplace = [];
    elseif 603==icycle,
        numchan = 19;
        for ii = 1:numchan,        nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );        fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\AhmetOmurtag0509\NIRS\';
        fpath2 = 'Finger\';
        fname = 'NIRS-2014-05-09_002';
        fpatheeg = 'F:\LabData\LabData\AhmetOmurtag0509\EEG\';
        fnameeeg = 'aomurtagfinger_20140509_173559';
        badchans = {'F8' 'T8' 'P8'};
        nirschanreplace = [];
    elseif 703==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\SarpOzture0527\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-05-27_001';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0527\EEG\';
        fnameeeg = 'TASK\sarptask_20140527_190028';
        badchans = { 'T8' };
        nirschanreplace = [];
        %{
        T8: sharp peaks in time series wl1,2. After convert/filter looks ok
        but PSD not good
        %}
    elseif 704==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\SarpOzture0527\NIRS\';
        fpath2 = 'Resting State\';
        fname = 'NIRS-2014-05-27_002';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0527\EEG\';
        fnameeeg = 'Resting State\sarpresting_20140527_191811';
        badchans = {};
        nirschanreplace = [];
    elseif 803==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\SarpOzture0601\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-06-01_001';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0601\EEG\';
        fnameeeg = 'sarptasksa_20140601_183535';
        badchans = {'T6' 'T8' };
        nirschanreplace = [];
        
        %{
        T6: serious jumps in time series wl1,wl2. after conversion and
        filter it's seriously messed up.
        T8: sharp peaks in time series wl1,2. After convert/filter looks ok
        but PSD not good
        FP1: jumps in time series wl1,2. After convert/filter looks ok. PSD
        ok.
        %}
        
    elseif 804==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\SarpOzture0601\NIRS\';
        fpath2 = 'Resting State\';
        fname = 'NIRS-2014-06-01_003';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0601\EEG\';
        fnameeeg = 'sarprestingsafinal_20140601_185422';
        badchans = { };
        nirschanreplace = [];
        %{
        T8: sharp peaks in time series wl1,2. After convert/filter looks
        ok. psd ok.
        %}
    elseif 903==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\Tom0604\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-06-04_002';
        fpatheeg = 'F:\LabData\LabData\Tom0604\EEG\';
        fnameeeg = 'tomtaskFINAL_20140604_161045';
        badchans = { };
        nirschanreplace = [];
        %{
        all chans have unusual round PSD peak at 7.5Hz
        %}
        
    elseif 904==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );    end
        fpath = 'F:\LabData\LabData\Tom0604\NIRS\';
        fpath2 = 'Resting\';
        fname = 'NIRS-2014-06-04_003';
        fpatheeg = 'F:\LabData\LabData\Tom0604\EEG\';
        fnameeeg = 'tomtRESTING_20140604_163208';
        badchans = { };
        nirschanreplace = [];
        %{        
        all chans have unusual round PSD peak at 7.5Hz
        all chans look artifacty first 2 mins
        %}
        
    elseif 1003==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            % fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\Lu0609\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-06-09_002';
        fpatheeg = 'F:\LabData\LabData\Lu0609\EEG\';
        fnameeeg = 'lufinaltask_20140609_173506';
        badchans = { };
        nirschanreplace = [];
        %{
    case 10,chanlabel = 'CZ';
    case 12,chanlabel = 'C4';
    case 13,chanlabel = 'O2';
        %}    
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
        badchans = { };
        nirschanreplace = [];
        %{        
    case 10,chanlabel = 'CZ';
    case 12,chanlabel = 'C4';
    case 13,chanlabel = 'O2';
        %}    
    end  
    %% GET / PREPROCESS EEG DATA
    % PREPROCESS params
    % data clip
    dtbeg = 10; % sec
    dtend = 10; % sec
    file1 = [fpatheeg fnameeeg '.edf'];
    [eegchanlabels, eegchanlocs, xxeegf, evteeg, Fs] = getpreprocess_eeg( file1, dtbeg, dtend, filterlofreq, filterhifreq, filternotch );
    xxeegf = -xxeegf; 
    if 0
        for ii = 1:numchan,
            ibeg = 1; iend = 10.*Fs;            chan2 = upper ( eegchanlabels(ii) );
            figure(980); clf; plot(  xxeegf(ibeg:iend,ii) ); title(['EEG  ' chan2]);
            figure(982); clf; nw = 3.5; pmtm( xxeegf(:, ii),  nw, [], Fs); title(chan2);
            pause(1);
        end
    end
    %% GET / PREPROCESS NIRS DATA
    % Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
    filterlofreq = .01; % Hz
    filterhifreq = .08; % Hz
    filternotch = 0; % 0: no notch;[
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
    
    % REPLACEMENT CORRECTION:
    [nrows,ncols] = size(nirschanreplace);
    if 0<nrows,
        for irow = 1:nrows,
            chanindxtobereplaced = nirschanreplace(irow,1);
            chanindxthatreplaces = nirschanreplace(irow,2);
            xx_wl1(:, chanindxtobereplaced) = xx_wl1(:, chanindxthatreplaces);
            xx_wl2(:, chanindxtobereplaced) = xx_wl2(:, chanindxthatreplaces);
        end
    end
    
    iflagdat = 13;
    %{
1: raw detector data,
2: filtered detector data
11: raw detector data converted to concentr GOOD
12: filtered detector data converted to concentr LEAST GOOD
13: raw detector data converted to concentr then filtered  BEST GOOD
    %}
    fignum=fignum+iflagdat;
    
    if 1==iflagdat,
        xxnirsf = xx_wl1;
        xxnirsfR = xx_wl2;
    elseif 2==iflagdat,
        xxnirsf = mybutterandiirnotchfilters( xx_wl1, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
        xxnirsfR = mybutterandiirnotchfilters( xx_wl2, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
    elseif 11==iflagdat,
        [xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
        xxnirsf = xxo;
        xxnirsfR = xxr;
    elseif 12==iflagdat,
        xx_wl1f = mybutterandiirnotchfilters( xx_wl1, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
        xx_wl2f = mybutterandiirnotchfilters( xx_wl2, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
        [xxo, xxr] = mybeerlambert2( xx_wl1f, xx_wl2f );
        xxnirsf = xxo;
        xxnirsfR = xxr;
    elseif 13==iflagdat,
        [xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
        xxnirsf = mybutterandiirnotchfilters( xxo, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
        xxnirsfR = mybutterandiirnotchfilters( xxr, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
    end;
    %% VIEW NIRS DATA
    if 1
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
            logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.3 );
            plot(F,logPxx); axis([0 .3 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
                        
            figure(470);clf;
            nw = 3.5; [Pxx,F] = pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
            logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.3 );
            plot(F,logPxx); axis([0 .3 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
               
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
    
    %figure(330); clf; plot(tteegclip, xxeegfclip(:,1), ttnirsclip, -50 + xxnirsfclip(:, 1) )
    %figure(332); clf;  plot(tteegclip, xxeegfclip(:,1), ttnirsclip, -50 + xxnirsfclip(:, 1), [ttnirsclip(nnn) ttnirsclip(nnn)], [-50 50],'r*' )
    %% TMP for kkk=1:19,plot(tteegclip, xxeegfclip(:,kkk)); title(num2str(kkk)); axis([0 30 -200 200]);pause(1);end
    
    %% save data
    xxeegfclip0 = xxeegfclip;
    tteegclip0 = tteegclip;
    xxnirsfclip0 = xxnirsfclip;
    xxnirsfRclip0 = xxnirsfRclip;
    ttnirsclip0 = ttnirsclip;
    %% EEG POWER            
    % calc power time course
    num = nweegpower-1; denom = nweegpower; noverlap = floor( nweegpower.* (num./denom) );
    freqs = [];
    nbegeegpower = floor(nweegpower./2);
    nendeegpower = length(xxeegfclip) - floor(nweegpower./2);
    xxeegpower = zeros(nendeegpower-nbegeegpower+1,numchan);
    
    Powarray = cell(numchan,1);
    
    for kchan = 1:numchan,
        fprintf('Calculating eeg power, chan %d/%d\n',kchan,numchan);
        [SS,Freq,Tim,Pow] = spectrogram(xxeegfclip(:,kchan),nweegpower,noverlap,freqs,Fs);
        Powarray{kchan} = Pow;
    end
    
    for ifreqs = [1:6],
                
        % start with fresh data
        xxeegfclip = xxeegfclip0;
        tteegclip = tteegclip0;
        xxnirsfclip = xxnirsfclip0;
        xxnirsfRclip = xxnirsfRclip0;
        ttnirsclip = ttnirsclip0;
        
        if 1==ifreqs,
            freqlolimit = 0; freqhilimit = 70;
        elseif 2==ifreqs,
            freqlolimit = 0; freqhilimit = 4;
        elseif 3==ifreqs,
            freqlolimit = 4; freqhilimit = 8;
        elseif 4==ifreqs,
            freqlolimit = 8; freqhilimit = 12;
        elseif 5==ifreqs,
            freqlolimit = 12; freqhilimit = 30;
        elseif 6==ifreqs,
            freqlolimit = 30; freqhilimit = 70;
        end
        
        for kchan = 1:numchan,
            fprintf('Assigning eeg power, chan %d/%d freqs %f %f\n',kchan,numchan,freqlolimit,freqhilimit);
            indxs = find( freqlolimit<=Freq & Freq<=freqhilimit );
            Pow = Powarray{kchan};
            %Pow = 10.*log10(Pow); % converting to dB has little effect on
            %lagged corr
            pa = sum( Pow(indxs,:), 1 );
            xxeegpower(:,kchan) = pa;
        end
        %tb=120; te=130; figure(7790);clf; plot(Tim,xxeegpower(:,kchan));axis([tb te 0 120])
        %figure(7110);clf; zz = 10.*log10(Pow); maxz = max(max(zz)); minz = min(min(zz)); surf(Tim,Freq,zz,'edgecolor','none'); axis([tb te 0 60 minz maxz] );view(0,90);xlabel('Time (s)'); ylabel('Hz');
        
        % assign eeg power time course to the eeg data
        xxeegfclipclip = xxeegpower;
        
        % cut eeg time to match eeg power
        tteegclipclip =  tteegclip( 1 : length(xxeegpower) );
        
        % cut nirs data
        indxbegnirs = ceil(nbegeegpower .* Fsnirs./Fs);
        indxendnirs = length(xxnirsfclip) - indxbegnirs + 1;
        xxnirsfclipclip = xxnirsfclip( indxbegnirs:indxendnirs, : );
        xxnirsfRclipclip = xxnirsfRclip( indxbegnirs:indxendnirs, : );

        % cut nirs time
        ttnirsclipclip = ttnirsclip(  indxbegnirs : indxendnirs );
        
        % assign clipped as new data
        xxeegfclip = xxeegfclipclip;
        tteegclip = tteegclipclip;
        xxnirsfclip = xxnirsfclipclip;
        xxnirsfRclip = xxnirsfRclipclip;
        ttnirsclip = ttnirsclipclip;
        
        clear xxeegfclipclip tteegfclipclip xxnirsfclipclip ttnirsfclipclip
        
        %% DOWNSAMPLE (INTERPOLATE) EEG TO NIRS DATA POINTS
        xxeegfclipds = zeros( length(ttnirsclip), numchan );
        for ii = 1:numchan,
%            xxeegfclipds(:,ii) = interp1( tteegclip, xxeegfclip(:,ii), ttnirsclip,'v5cubic' );            
            FF = griddedInterpolant( tteegclip, xxeegfclip(:,ii), 'spline','none' );
            xxeegfclipds(:,ii) = FF(ttnirsclip);
        end
        %% VIEW EEG
        if 0
            for ii = 1:numchan,
                ibeg = 1; iend = 10.*Fs;            chan2 = upper ( eegchanlabels(ii) );
                figure(980); clf; plot(  xxeegfclipds(:,ii) ); title(['EEG Power ' chan2]);
                figure(982); clf; pmtm( xxeegfclipds(:,ii ) ); title(chan2);
                pause(1);
            end
        end
        %%
        for ii = 1:numchan,
            if find(1==isnan(xxnirsfclip(:,ii))) ,
                warning('isnan 16055');
            end
        end
        %% LAGGED CORR
        Nwmax = 100000000; % max num windows to use (large number implies use max available)
        %nlags = ceil( Fsnirs.*(Tw-1) ); % max number of lags to compute (samples)
        nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
        
        numsampW = floor( Fsnirs.*Tw );
        nstds = [];
        
        maxdim = min( Nwmax, floor( length(xxnirsfclip)./numsampW ) );
        
        % first window
        iw = 1;
        ibeg = (iw-1).*numsampW + 1;
        iend = ibeg + numsampW - 1;
        
        XCFmean = []; XCFRmean = []; XCF3mean = []; XCFmean_surr = []; XCFmean_surr2 = [];
        XCFmean_shuff = [];
        XCFmean_delayed = [];
        
        fprintf( 'Calculating lagged correlation maxima for all windows all channel pairs\n' );
        
        % loop through windows
        while iw<=Nwmax & iend<=length(xxnirsfclip),
            
            dat1 = xxeegfclipds( ibeg:iend, : );
            dat2 = xxnirsfclip( ibeg:iend, : );
            dat2R = xxnirsfRclip( ibeg:iend, : );
            
            % delayed
            ndelay = ceil( 5 .*Fsnirs );
            dat1_delayed = dat1; dat1_delayed(ndelay+1:end,:) = dat1(1:end-ndelay,:);
            
            % surrogate data (shuffled time points)
            dat2_shuff = zeros( size(dat2) );
            for kchan=1:numchan, dat2_shuff(:,kchan) = dat2(randperm(length(dat2)),kchan); end
            
            % surrogate data obtained from dat1 with known delayed corr
            % gaussian
            winsiz = 60; tdely1 = 5; tdely2 = 5;
            ttt = [1:ceil(winsiz.*Fsnirs)]./Fsnirs;
            mu = ceil(winsiz./2+tdely1); ss = ceil(4); fff1 = exp( -(ttt-mu).^2./(2.*ss) );
            mu = ceil(winsiz./2+tdely2); ss = ceil(2); fff2 = -exp( -(ttt-mu).^2./(2.*ss) );
            % gamma function
            midsz = ceil(1..*length(ttt)./2); tttmidsz = fliplr( ttt(1:midsz) );
            fff3 = zeros( size(ttt) );
            alph = 3; fff3(1:midsz) = fff3(1:midsz) + (tttmidsz).^(alph-1).*exp(-(tttmidsz./8));
            fff4 = zeros( size(ttt) );
            alph = 3; fff4(1:midsz) = -(tttmidsz).^(alph-1).*exp(-(tttmidsz./8));
            
            fff1 = fff1./abs(sum(fff1));
            fff2 = fff2./abs(sum(fff2));
            fff3 = fff3./abs(sum(fff3));
            fff4 = fff4./abs(sum(fff4));
            
            %fff3 = fff1+ fff3; fff4 = fff2+ fff4;
            fff3 = fff1; fff4 = fff2;
            
            if 1==iw, figure(2222);clf;plot(ttt,fff1,ttt,fff2,ttt,fff3,ttt,fff4); end
            
            dat_surr = zeros( size(dat2) );
            dat_surr2 = zeros( size(dat2) );
            for kchan=1:numchan, dat_surr(:,kchan) = conv(dat1(:,kchan),fff3,'same'); end
            for kchan=1:numchan, dat_surr2(:,kchan) = conv(dat1(:,kchan),fff4,'same'); end
            
            if 1
                dat1 = mynormalize( dat1 );
                dat2 = mynormalize( dat2 );
                dat2R = mynormalize( dat2R );
                dat2_shuff = mynormalize( dat2_shuff );
                dat_surr = mynormalize( dat_surr );
                dat_surr2 = mynormalize( dat_surr2 );
                dat1_delayed = mynormalize( dat1_delayed);
            end
            
            for ii = 1:numchan, % NIRS
                chan1 = upper( nirschanlabels(ii) );
                
                indxbadchan = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
                if 0~=indxbadchan, continue; end
                
                %if ~strcmp( chan1, 'C3' ) & ~strcmp( chan1, 'C4' )  , continue; end
                %if strcmp( chan1, 'FP1' ) | strcmp( chan1, 'FP2' )  , continue; end
                if flag_useprefrontalonly&~isprefrontal( char(chan1) ), continue; end
                
                % Tinh good channels
                %if ~(1==ii | 2==ii | 3==ii | 4==ii | 9==ii | 16==ii | 18==ii | 19==ii), continue; end
                
                for jj = 1:numchan, % EEG
                    chan2 = upper ( eegchanlabels(jj) );
                    
                    indxbadchan = findstrincellarray( badchans, upper(chan2) ); % return zero if chan1 is NOT a bad channel
                    if 0~=indxbadchan, continue; end
                    
                    if strcmp( chan1, chan2 ),
                        fprintf( 'iw %5d/%5d  %6s %6s \n', iw, maxdim, char(chan1), char(chan2) );
                        
                        % varlags = [1:floor(nlags./2)]; vars = calc_crosscorrvariance( dat2(:,jj), dat1(:,ii), nlags, varlags );
                        
                        if find(1==isnan(dat2(:,ii))) | find(1==isnan(dat1(:,jj))),
                            warning('isnan 2010');
                            adsffsda=32425;
                        end
                        
                        [ XCF, Lags ] = xcorr(  dat2(:,ii), dat1(:,jj),  nlags, 'coeff' );
                        %[ XCF ] = laggedmutualinfo(  dat2(:,ii)', dat1(:,jj)',  nlags)';
                        [ XCFR, Lags ] = xcorr(   dat2R(:,ii), dat1(:,jj), nlags, 'coeff' );
                        %[ XCFR ] = laggedmutualinfo(  dat2R(:,ii)', dat1(:,jj)',  nlags)';
                        [ XCF3, Lags ] = xcorr( dat2R(:,ii),  dat2(:,ii),  nlags, 'coeff' );
                        [ XCFshuff, Lags ] = xcorr(  dat2_shuff(:,jj), dat1(:,ii),  nlags, 'coeff' );
                        [ XCFs, Lags ] = xcorr( dat_surr(:,jj),  dat1(:,ii),  nlags, 'coeff' );
                        [ XCFs2, Lags ] = xcorr( dat_surr2(:,jj),  dat1(:,ii),  nlags, 'coeff' );
                        [ XCFd, Lags ] = xcorr(  dat1_delayed(:,jj),  dat1(:,ii), nlags, 'coeff' );
                        
                        
                        XCFmean = [ XCFmean, XCF ];
                        XCFRmean = [ XCFRmean, XCFR ];
                        XCF3mean = [ XCF3mean, XCF3 ];
                        XCFmean_shuff = [ XCFmean_shuff, XCFshuff];
                        XCFmean_surr = [ XCFmean_surr, XCFs ];
                        XCFmean_surr2 = [ XCFmean_surr2, XCFs2 ];
                        XCFmean_delayed = [ XCFmean_delayed, XCFd ];
                        
                        if false, % plot correlation vs. lag
                            figure(880); clf;
                            plot( Lags./Fsnirs, XCF, 'r' ); hold on;
                            plot( Lags./Fsnirs, XCFR, 'g' ); hold on;
                            plot( varlags./Fsnirs, vars, 'k' );
                            maxlim = max(Lags./Fsnirs);
                            axis([-maxlim  maxlim  -.5  .5]);
                            grid on;
                            titlestr = [  'iw ' num2str(iw) ', chans ' char(chan1) ' ' char(chan2) ];
                            title( titlestr );
                            pause(1);
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
        
        
        %%
        XCFmean = mean( XCFmean, 2 );
        XCFRmean = mean( XCFRmean, 2 ); XCF3mean = mean( XCF3mean, 2 );
        XCFmean_surr = mean( XCFmean_surr, 2 ); XCFmean_surr2 = mean( XCFmean_surr2, 2 );
        XCFmean_shuff = mean( XCFmean_shuff, 2);
        XCFmean_delayed = mean( XCFmean_delayed, 2);
        %% PSD
        if 0
            for ii = 1:numchan,
                chan1 = upper( nirschanlabels(ii) );
                if ~isprefrontal( char(chan1) ), continue; end
                for jj = 1:numchan,
                    chan2 = upper ( eegchanlabels(jj) );
                    if strcmp( chan1, chan2 ),
                        
                        dat1 = xxeegfclipds( : , : );
                        dat2 = xxnirsfclip( : , : );
                        
                        dat1 = mynormalize( dat1 );
                        dat2 = mynormalize( dat2 );
                    end
                end
            end
        end
        %% PLOT
        gry = [.9 .9 .9]; %gry = [.8 .9 .8];
        fh = figure(fignum+icycle); clf;
%        set(fh,'position',[10 130 280 180]); set(fh, 'color','w' ); 
        set(fh,'position',[10 130 840 540]); set(fh, 'color','w' );

        % smooth
        ns = 5;
        XCFmeans = mysmooth(XCFmean,ns);
        XCFRmeans = mysmooth(XCFRmean,ns);
        XCF3means = mysmooth(XCF3mean,ns);
        XCFmeans_surr  = mysmooth(XCFmean_surr,ns);
        XCFmeans_surr2  = mysmooth(XCFmean_surr2,ns);
        XCFmeans_shuff = mysmooth(XCFmean_shuff,ns);
        %pp1 = plot( Lags./Fsnirs, .7.*XCFmeans_surr,'color', gry ,'linestyle','--','linewidth', linsz2); hold on;
        %pp1 = plot( Lags./Fsnirs, XCFmeans_surr2,'color', 'b' ,'linestyle','--','linewidth', linsz2); hold on;
        pp = plot( Lags./Fsnirs, XCFmeans, 'color', 'r' , 'linewidth', linsz); hold on;
        pp = plot( Lags./Fsnirs, XCFRmeans,'color', 'b' , 'linewidth', linsz); hold on;
        %pp = plot( Lags./Fsnirs, XCF3means,'color', 'k' , 'linewidth', 2); hold on;
        %pp = plot( Lags./Fsnirs, XCFmeans_shuff,'color', 'k' , 'linewidth', .1); hold on;
        axis(  [xminlim  xmaxlim  yminlim  ymaxlim]  );
        %set( gca, 'fontsize', fontsztic );
        if 0
            if flagxtics, xlabel('Time Lag (s)', 'fontsize', fontsz); set(gca,'xTick',[-25:10:25],'fontsize',fontsz);
            else, set(gca,'xTick',[]); end
            if flagytics, ylabel('Correlation', 'fontsize', fontsz); set(gca,'yTick',[yminlim:.01:ymaxlim],'fontsize',fontsz);
            else, set(gca,'yTick',[]); end
        end
        grid on;
        mystrn = [ '_win_' num2str( Tw ) '_freqs' num2str( freqlolimit ) '-' num2str( freqhilimit ) '_icyc_' num2str(icycle)];
        titlestrn = {...
            [fpath fpath2];...
            [fname ];...
            mystrn;...
            };
        
        if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b','fontsize',9); end
        drawnow;
        %%
        if flag_print_fig,
            set( gcf, 'PaperPositionMode', 'auto')
            print(  '-dtiff', '-zbuffer', '-r600',  [figspath  figname mystrn '.tiff'] );
        end
        %% PLOT
        if 0
            fontsz = 22; fontsztic = 16; linsz = 3;
            fh = figure(2000+fignum+icycle); clf;
            set(fh, 'color','w' );
            % smooth
            ns = 5;
            XCFmeans_surr  = mysmooth(XCFmean_surr,ns);
            XCFmeans_surr2  = mysmooth(XCFmean_surr2,ns);
            XCFmeans_delayed  = mysmooth(XCFmean_delayed,ns);
            
            %    xminlim = -Tw; xmaxlim = Tw; yminlim = -.04; ymaxlim = .04;
            %xminlim = -5; xmaxlim = 5; yminlim = -1.05; ymaxlim = 1.05;
            pp3 = plot( Lags./Fsnirs, XCFmeans_delayed,'color', 'g' , 'linewidth', .1); hold on;
            pp1 = plot( Lags./Fsnirs, XCFmeans_surr,'color', 'r' , 'linewidth', .1); hold on;
            pp2 = plot( Lags./Fsnirs, XCFmeans_surr2,'color', 'b' , 'linewidth', .1); hold on;
            axis(  [xminlim  xmaxlim  yminlim  ymaxlim]  );
            %   xlabel('Time Lag (s)', 'fontsize', fontsz); ylabel('Correlation', 'fontsize', fontsz);
            set( gca, 'fontsize', fontsztic );
            %    set(gca,'xTick',[-15:5:15]);    set(gca,'yTick',[yminlim:ymaxlim./2:ymaxlim]);
            
            titlestrn = {...
                [fpath fpath2];...
                [fname ];...
                [ 'Window size (s): ' num2str( Tw ) ];...
                [ '_freqs' num2str( freqlolimit ) '-' num2str( freqhilimit ) '_icyc_' num2str(icycle)];...
                };
            
            if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b','fontsize',9); end
            grid on
            
        end
        %% PLOT
        gry = [.9 .9 .9]; %gry = [.8 .9 .8];
        fh = figure(fignum+icycle+1); clf;
        set(fh,'position',[10 555 475 350]); set(fh, 'color','w' );
        % smooth
        ns = 5;
        XCFmeans = mysmooth(XCFmean,ns);
        XCFRmeans = mysmooth(XCFRmean,ns);
        XCF3means = mysmooth(XCF3mean,ns);
        XCFmeans_surr  = mysmooth(XCFmean_surr,ns);
        XCFmeans_surr2  = mysmooth(XCFmean_surr2,ns);
        XCFmeans_shuff = mysmooth(XCFmean_shuff,ns);
        %pp1 = plot( Lags./Fsnirs, .7.*XCFmeans_surr,'color', gry ,'linestyle','--','linewidth', linsz2); hold on;
        %pp1 = plot( Lags./Fsnirs, XCFmeans_surr2,'color', 'b' ,'linestyle','--','linewidth', linsz2); hold on;
        %pp = plot( Lags./Fsnirs, XCFmeans, 'color', 'r' , 'linewidth', linsz); hold on;
        %pp = plot( Lags./Fsnirs, XCFRmeans,'color', 'b' , 'linewidth', linsz); hold on;
        pp = plot( Lags./Fsnirs, XCF3means,'color', 'k' , 'linewidth', 2); hold on;
        %pp = plot( Lags./Fsnirs, XCFmeans_shuff,'color', 'k' , 'linewidth', .1); hold on;
        axis(  [xminlim  xmaxlim  -.3  .3]  );
        %set( gca, 'fontsize', fontsztic );
        if 0
            if flagxtics, xlabel('Time Lag (s)', 'fontsize', fontsz); set(gca,'xTick',[-25:10:25],'fontsize',fontsz);
            else, set(gca,'xTick',[]); end
            if flagytics, ylabel('Correlation', 'fontsize', fontsz); set(gca,'yTick',[yminlim:.01:ymaxlim],'fontsize',fontsz);
            else, set(gca,'yTick',[]); end
        end
        grid on;
        titlestrn = {...
            [fpath fpath2];...
            [fname ];...
            [ 'Window size (s): ' num2str( Tw ) ];...
            };
        
        if idisplaytitle, title( titlestrn , 'Interpreter','none','color','b','fontsize',9); end
        
    end % EEG POWER
end % ICYCLE
