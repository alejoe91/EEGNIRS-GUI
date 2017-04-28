function datacells = geteegnirsdataALE( ...
    experimentid,...
    dtbegeeg,...
    dtendeeg,...
    filterlofreqeeg,...
    filterhifreqeeg,...
    filternotcheeg,...
    numsampleswindoweegpower,...
    filterlofreqnirs,...
    filterhifreqnirs ,...
    filternotchnirs,...
    flag_viewnirs, ...
    flag_vieweeg, ...
    tw, ...
    overlap...
    )

datacells = cell( length(experimentid),  22 );

%{
 DATA STRUCTURE:

PREPROCESSED
xxeegf
xxnirsf
xxnirsfR

SYNCED
xxeegfclip
tteegclip
xxnirsfclip
xxnirsfRclip
ttnirsclip

EEG resampled to NIRS time points
xxeegfclipds

EEG bandpower  ( time x chans x freqbands )
xxeegfclipbandpower

EEG bandpower resampled to NIRS time points ( time x chans x freqbands )
xxeegfclipbandpowerds

EVENT TIMES (data clipped at first event, thus first event time is 0)
eventtimes

%}

count = 1;

%folder0 = 'Z:\xys\';

for icycle = experimentid,
    
    if 703==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            %   fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\SarpOzture0527\NIRS\';
        
  %              fpath =  [ folder0 '\LabData\SarpOzture0527\NIRS\' ];

        fpath2 = 'Task\';
        fname = 'NIRS-2014-05-27_001';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0527\EEG\';
        fnameeeg = 'TASK\sarptask_20140527_190028';
        badchans = { 'T8' };
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
        %{
        T8: sharp peaks in time series wl1,2. After convert/filter looks ok
        but PSD not good
        %}
    elseif 704==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            %  fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\SarpOzture0527\NIRS\';
        fpath2 = 'Resting State\';
        fname = 'NIRS-2014-05-27_002';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0527\EEG\';
        fnameeeg = 'Resting State\sarpresting_20140527_191811';
        badchans = {};
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
    elseif 803==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            %   fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\SarpOzture0601\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-06-01_001';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0601\EEG\';
        fnameeeg = 'sarptasksa_20140601_183535';
        badchans = {'T6' 'T8' };
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
        
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
            % fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\SarpOzture0601\NIRS\';
        fpath2 = 'Resting State\';
        fname = 'NIRS-2014-06-01_003';
        fpatheeg = 'F:\LabData\LabData\SarpOzture0601\EEG\';
        fnameeeg = 'sarprestingsafinal_20140601_185422';
        badchans = { };
        numsrc = 24;
        numdet = 16;
        nirschanreplace = [];
        indxchanssubset = [];
        %{
        T8: sharp peaks in time series wl1,2. After convert/filter looks
        ok. psd ok.
        %}
    elseif 903==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            % fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\Tom0604\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-06-04_002';
        fpatheeg = 'F:\LabData\LabData\Tom0604\EEG\';
        fnameeeg = 'tomtaskFINAL_20140604_161045';
        badchans = { };
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
        %{
        all chans have unusual round PSD peak at 7.5Hz
        %}
        
    elseif 904==icycle,
        numchan = 16;
        for ii = 1:numchan,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair( ii, ii );
            %     fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\Tom0604\NIRS\';
        fpath2 = 'Resting\';
        fname = 'NIRS-2014-06-04_003';
        fpatheeg = 'F:\LabData\LabData\Tom0604\EEG\';
        fnameeeg = 'tomtRESTING_20140604_163208';
        badchans = { };
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
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
        badchans = { 'CZ' };
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
        %{
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
        badchans = {'CZ' };
        nirschanreplace = [];
        numsrc = 24;
        numdet = 16;
        indxchanssubset = [];
        %{
        %}
    elseif 1102==icycle,
        numchanmax = 24;
        for ii = 1:numchanmax,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair_multi( icycle, ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\Adada0626\fNIRS\';
        fpath2 = 'Resting State\';
        fname = 'NIRS-2014-06-26_001';
        fpatheeg = 'F:\LabData\LabData\Adada0626\EEG\';
        fnameeeg = 'Resting State\zxzxczxc_20140626_170538';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 32;
        numdet = 24;
        indxchanssubset = [1:12 17:24];
        %{
       NIRS
        heartrate xxo PSD peak not sharp:         T7 P7 CZ T8
        POOR CALIB: C3 P7 O1
        EEG
        artifact from NIRS visible in P3
        %}
    elseif 1103==icycle,
        numchanmax = 24;
        for ii = 1:numchanmax,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair_multi( icycle, ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\Adada0626\fNIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-06-26_002';
        fpatheeg = 'F:\LabData\LabData\Adada0626\EEG\';
        fnameeeg = 'Task\zxzxczxctask_20140626_172455';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 32;
        numdet = 24;
        indxchanssubset = [1:12 17:24];
        %{
       NIRS
        heartrate xxo PSD peak not sharp:          P7  T8
        POOR CALIB: C3 P7 O1
        EEG
        nirs artifact not visible
        %}
    elseif 1112==icycle,
        numchanmax = 24;
        for ii = 1:numchanmax,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair_multi( icycle, ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\HalehAghajani0702\NIRS\';
        fpath2 = 'Resting\';
        fname = 'NIRS-2014-07-02_001';
        fpatheeg = 'F:\LabData\LabData\HalehAghajani0702\EEG\';
        fnameeeg = 'hallrest_20140702_203131';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 32;
        numdet = 24;
        indxchanssubset = [1:12 17:24];
        %{
       NIRS        
        EEG
        %}
    elseif 1113==icycle,
        numchanmax = 24;
        for ii = 1:numchanmax,
            nirschanlabels{ii} = get_chanlabel_from_srcdetpair_multi( icycle, ii, ii );
            fprintf( '%d %d  %s\n', ii, ii, nirschanlabels{ii} );
        end
        fpath = 'F:\LabData\LabData\HalehAghajani0702\NIRS\';
        fpath2 = 'Task\';
        fname = 'NIRS-2014-07-02_003';
        fpatheeg = 'F:\LabData\LabData\HalehAghajani0702\EEG\';
        fnameeeg = 'halltask_20140702_204922';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 32;
        numdet = 24;
        indxchanssubset = [1:12 17:24];
        %{
       NIRS        
        EEG
        %}
    elseif 1122==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\AhmetOmurtag0722\PronounceVFT\';
        % nirs
        fpath = [ parentfold 'NIRS\' ];
        fpath2 = '';
        fname = 'NIRS-2014-07-22_002';

        fpatheeg = [ parentfold 'EEG\'];
        fnameeeg = 'ahmethocanewtaskfinal_20140722_163730';
        badchans = { 'CZ' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];
        %{
       NIRS        
        EEG - CZ is bad
        %}
    elseif 1123==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\AhmetOmurtag0722\WritingVFT\';
        % nirs
        fpath = [ parentfold 'fNIRS\' ];
        fpath2 = '';
        fname = 'NIRS-2014-07-22_003';

        fpatheeg = [ parentfold 'EEG\'];
        fnameeeg = 'ahmethocanewwritingfinalfinal_20140722_165150';
        badchans = { 'CZ' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];
        %{
       NIRS        
        EEG - CZ is bad
        %} 
    elseif 1132==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'C:\Users\alessio\Documents\Ale\UH\Thesis\NIRS_EEG\Data\raw_data\AhmetOmurtag0801\';
        % nirs
        fpath = [ parentfold 'fNIRS\RestingState\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-01_002';

        fpatheeg = [ parentfold 'EEG\Resting State\'];
        fnameeeg = 'Omurtag3EMFResting_20140801_160523';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];    
    elseif 1133==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\AhmetOmurtag0801\';
        % nirs
        fpath = [ parentfold 'fNIRS\TaskVFT\VFT - Pronounce\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-01_003';

        fpatheeg = [ parentfold 'EEG\TaskVFT\Pronounce\'];
        fnameeeg = 'Omurtag3EMFVFTpronounce_20140801_163220';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];      
    elseif 1134==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\AhmetOmurtag0801\';
        % nirs
        fpath = [ parentfold 'fNIRS\TaskVFT\VFT- Writing\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-01_004';

        fpatheeg = [ parentfold 'EEG\TaskVFT\Writing\'];
        fnameeeg = 'Omurtag3EMFVFTwriting_20140801_164130';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];        
    elseif 1135==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
  %findstr:returns a number if it finds S2 in S1; isempty returns 0 if the vector is empty --> it enters the if else if tmpstring!='UNUSED'
                if isempty(findstr(tmpstrng, 'UNUSED')), 
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'C:\Users\alessio\Documents\Ale\UH\Thesis\NIRS_EEG\Data\raw_data\AhmetOmurtag0801\';
        % nirs
        fpath = [ parentfold 'fNIRS\Artifacts\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-01_005';

        fpatheeg = [ parentfold 'EEG\Artifact\'];
        fnameeeg = 'Omurtag3EMFVFTartifacts_20140801_165430';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];           
    elseif 1142==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\SarpOzture0822\';
        % nirs
        fpath = [ parentfold 'NIRS\Resting State\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-22_001';

        fpatheeg = [ parentfold 'EEG\Resting State\'];
        fnameeeg = 'RestingStateSarp_20140822_184420';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];              
    elseif 1143==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\SarpOzture0822\';
        % nirs
        fpath = [ parentfold 'NIRS\VFT-Pronounce\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-22_004';

        fpatheeg = [ parentfold 'EEG\VFT-Pronounce\'];
        fnameeeg = 'TaskpronouncefinalfinalSarp_20140822_191057';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];                 
    elseif 1144==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\SarpOzture0822\';
        % nirs
        fpath = [ parentfold 'NIRS\VFT-Writing\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-22_005';

        fpatheeg = [ parentfold 'EEG\VFT-Writing\'];
        fnameeeg = 'TaskwritingfinalSarp_20140822_192059';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];                   
    elseif 1145==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\SarpOzture0822\';
        % nirs
        fpath = [ parentfold 'NIRS\Artifacts\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-22_006';

        fpatheeg = [ parentfold 'EEG\Artifacts\'];
        fnameeeg = 'ArtifactSarp_20140822_193044';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];             
    elseif 1152==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\Tom0825\';
        % nirs
        fpath = [ parentfold 'fNIRS\RestingState\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-25_001';

        fpatheeg = [ parentfold 'EEG\RestingState\'];
        fnameeeg = 'tomrestingstate_20140825_165835';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];                
    elseif 1153==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\Tom0825\';
        % nirs
        fpath = [ parentfold 'fNIRS\VFT-Pronounce\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-25_002';

        fpatheeg = [ parentfold 'EEG\VFT-Pronounce\'];
        fnameeeg = 'tomVFTpronounce_20140825_171757';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];                 
    elseif 1154==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\Tom0825\';
        % nirs
        fpath = [ parentfold 'fNIRS\VFT-Writing\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-25_003';

        fpatheeg = [ parentfold 'EEG\VFT-Writing\'];
        fnameeeg = 'tomVFTwriting_20140825_172650';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];                    
    elseif 1155==icycle,
        numchanmax = 24;
        ichancount = 1;
        for ii = 1:numchanmax,
            for jj = 1:numchanmax,
                tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, ii, jj );
                if isempty(findstr(tmpstrng, 'UNUSED')),
                    nirschanlabels{ichancount} = tmpstrng;
                    fprintf( 'chan %d,    src,det %d %d,    label %s\n', ichancount, ii, jj, nirschanlabels{ichancount} );
                    ichancount = ichancount + 1;
                end
            end
        end
        numchan = ichancount - 1;
        parentfold = 'F:\LabData\LabData\Tom0825\';
        % nirs
        fpath = [ parentfold 'fNIRS\Artifacts\' ];
        fpath2 = '';
        fname = 'NIRS-2014-08-25_004';

        fpatheeg = [ parentfold 'EEG\Artifacts\'];
        fnameeeg = 'tomVFTartifact_20140825_173525';
        badchans = { '' };
        nirschanreplace = [];
        numsrc = 20;
        numdet = 24;
        indxchanssubset = [1:numchan];    
    else,
        error('not supported')
    end
    
    fprintf('[%s]\n',[fpath fpath2 fname]  );
    fprintf('[%s]\n', [fpatheeg fnameeeg] );
    %% get EEG
    dtbeg = dtbegeeg;
    dtend = dtendeeg;
    % EEG  Bandpas filter and notch (normal bandpass range: 0.5-70 Hz)
    filterlofreq = filterlofreqeeg;
    filterhifreq = filterhifreqeeg; %
    filternotch = filternotcheeg;
    file1 = [fpatheeg fnameeeg '.edf'];
    [eegchanlabels, eegchanlocs, xxeegf, evteeg, Fs] = getpreprocess_eeg( file1, dtbeg, dtend, filterlofreq, filterhifreq, filternotch );
    numchaneeg = length(eegchanlabels);
    %% get NIRS
    %% GET / PREPROCESS NIRS DATA
    % Bandpass filter and notch (normal bandpass range: 0.01-.2 Hz)
    filterlofreq = filterlofreqnirs;
    filterhifreq = filterhifreqnirs; %
    filternotch = filternotchnirs;
    % GET sample rate
    Fsnirs = str2num( getfieldval( [fpath fpath2 fname '_config.txt'], 'SamplingRate', '=') );
    %  EVENTS filename
    filenirsevt = [ fpath  fpath2  fname  '.evt' ];
    
    % LOAD DATA
    mydat = load( [ fpath  fpath2  fname '.wl1' ]  );
    xx = fnnirs_reduce_data_3( icycle, numsrc, numdet, mydat );
    xx_wl1 = xx;
    mydat = load( [ fpath  fpath2  fname '.wl2' ]  );
    xx = fnnirs_reduce_data_3( icycle, numsrc, numdet, mydat );
    xx_wl2 = xx;
    % wl1--> red  wl2 --> ir
    
    if ~isempty( indxchanssubset ),
        xx_wl1 = xx_wl1(:, indxchanssubset );
        xx_wl2 = xx_wl2(:, indxchanssubset );
        nirschanlabels = nirschanlabels( indxchanssubset  );
        numchan = length( indxchanssubset );
    end
    
    % convert to Hb
    [xxo, xxr] = mybeerlambert2( xx_wl1, xx_wl2 );
    % preprocess
    xxnirsf = mybutterandiirnotchfilters( xxo, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
    xxnirsfR = mybutterandiirnotchfilters( xxr, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
    %% VIEW NIRS DATA
    if flag_viewnirs
        for ii = 1:numchan,
            chan1 =  upper( nirschanlabels{ii} );
            
            figure(330);clf;
            tt = [1:length(xxnirsf)]./Fsnirs;
            subplot(2,1,1),plot( tt, xxnirsf(:, ii)); title( [num2str(ii) ' ' chan1 ' preprocessed' ] , 'fontsize', 22);
            subplot(2,1,2),plot( tt, xxnirsfR(:, ii) );
            
            figure(332);clf;
            tt = [1:length(xxnirsf)]./Fsnirs;
            subplot(2,1,1),plot( tt, xx_wl1(:, ii)); title( [num2str(ii) ' ' chan1 ' wl1 wl2'] , 'fontsize', 22);
            subplot(2,1,2),plot( tt, xx_wl2(:, ii) );
            
            figure(440);clf;
            nw = 3.5; pmtm( xxnirsf(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' chan1  ' HBO preprocessed' ] , 'fontsize', 22);
            
            figure(450);clf;
            nw = 3.5; pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' chan1  ' HBR preprocessed' ] , 'fontsize', 22);
            
            figure(442);clf;
            nw = 3.5; pmtm( xxo(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' chan1  ' HBO raw' ] , 'fontsize', 22);
            
            figure(452);clf;
            nw = 3.5; pmtm( xxr(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' chan1 ' HBR raw' ] , 'fontsize', 22);
            
            figure(460);clf;
            nw = 3.5; [Pxx,F] = pmtm( xxnirsf(:, ii),  nw, [], Fsnirs);
            logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.8 );
            plot(F,logPxx); axis([0 .8 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
            title( [num2str(ii) '  -  ' chan1  ' HBO preprocessed' ] , 'fontsize', 22);
            
            figure(470);clf;
            nw = 3.5; [Pxx,F] = pmtm( xxnirsfR(:, ii),  nw, [], Fsnirs);
            logPxx = 10.*log10(Pxx); indxs = find( 0<F & F<.8 );
            plot(F,logPxx); axis([0 .8 min(logPxx(indxs)) max(logPxx(indxs))]);grid on;
            title( [num2str(ii) '  -  ' chan1 ' HBR preprocessed' ] , 'fontsize', 22);
            
            figure(494);clf;
            nw = 3.5; pmtm( xx_wl1(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' chan1  ' wl1 raw' ] , 'fontsize', 22);
            
            figure(496);clf;
            nw = 3.5; pmtm( xx_wl2(:, ii),  nw, [], Fsnirs);
            title( [num2str(ii) '  -  ' chan1  ' wl2 raw' ] , 'fontsize', 22);
            
            pause(.1);
            result = input('Press <Enter> to continue...');
        end
    end
    %% SYNC
    % EEG event sample numbers
    fprintf('Clipping data to sync\n');
    Vevthresh = 1000;
    %finds the samples above threshold (from 106501 to 21 points)
    nevteeg = remove_contiguous(  find( Vevthresh<evteeg )  );
    tteeg = [1:length(evteeg)]./Fs;
    tevteeg = nevteeg./Fs;
    %figure(110);clf; plot( tteeg, evteeg ); hold on; plot( tevteeg, ones(length(nevteeg)).*Vevthresh,'r.'); %axis([60 65 -500 1200])
    %Same machine sampling eeg-nirs-->only the first event needed for sync
    nbegeeg = nevteeg(1);
    
    % NIRS event sample numbers
    myevtdat = load( [ filenirsevt ]  );
    nevtnirs = myevtdat(:,1);
    ttnirs = [1:length(xxnirsf)]./Fsnirs;
    tevtnirs = nevtnirs./Fsnirs;
    %figure(210);clf; plot( ttnirs, xxnirsf(:,1) ); hold on; plot( tevtnirs, ones(length(nevtnirs)).*.04,'r.')
    nbegnirs = nevtnirs(1);
    
    % slippage ??
    tslip = tevteeg(1)-tevtnirs(1);
    %tevteeg
    %tevtnirs + tslip
    
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
    %xxeegfclipds = zeros( length(ttnirsclip), numchan );
    %for ii = 1:numchaneeg,
    %    FF = griddedInterpolant( tteegclip, xxeegfclip(:,ii), 'spline','none' );
    %    xxeegfclipds(:,ii) = FF(ttnirsclip);
    %end
    %% view eeg
    if flag_vieweeg
        for ii = 1:numchaneeg,
            chan1 = upper( eegchanlabels{ii} );
            nw=3.5;
            figure(2330);clf; plot( tteegclip, xxeegfclip(:,ii) );title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
            figure(2340);clf; pmtm( xxeegfclip(:,ii),nw,[],Fs);title( [num2str(ii) ' ' chan1] , 'fontsize', 22);
            pause(.1);
            result = input('Press <Enter> to continue...');
        end
    end
    %% preparations for taking out eeg modes
    if 0
        for ii = 1:numchaneeg,
            
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
    %nweegpower = numsampleswindoweegpower;
    nweegpower = floor(tw * Fs);
    %    nweegpower = 50; % (must be even) window size in num samples
    %num = nweegpower-1; denom = nweegpower;
    %why? nooverlap = num is not enough? to be reviewed
    %if (mod(nweegpower,2)>=1)
        noverlap = floor( (overlap / 100) * nweegpower);
    %else
    %    noverlap = floor( nweegpower - nweegpower/2);
    %end
    
    freqs = nweegpower;
    
    numfreqbands = 8;
    xxeegfclipbandpower = zeros( floor(length(tteegclip)/nweegpower), numchaneeg, numfreqbands );
    %xxeegfclipbandpowerds = zeros( length(ttnirsclip), numchaneeg, numfreqbands );
    
    for kchan = 1:numchaneeg,
        fprintf('Calculating eeg power, chan %d/%d\n',kchan,numchaneeg);
        %freqs = [];
        [SS,Freq,Tim,Pow] = spectrogram(xxeegfclip(:,kchan),nweegpower,noverlap,freqs,Fs);
        %Pow = 10.*log10(Pow); % converting to dB has little effect on lagged corr
        
        % Pow is shorter than the original time series because of window
        % size. Replicate first/last entry to pad the missing initial/final
        % segments
        %avecinitial = Pow(:,1); avecfinal = Pow(:,end);
        %Pow = [ repmat(avecinitial,1,nweegpower./2-1) Pow repmat(avecfinal,1,nweegpower./2) ];
        
        % loop over frequency bands
        for ifreqs = [1:numfreqbands],
            if 1==ifreqs, freqlolimit = 0; freqhilimit = 80;
            elseif 2==ifreqs, freqlolimit = 0; freqhilimit = 4;
            elseif 3==ifreqs, freqlolimit = 4; freqhilimit = 8;
            elseif 4==ifreqs, freqlolimit = 8; freqhilimit = 12;
            elseif 5==ifreqs, freqlolimit = 12; freqhilimit = 30;
            elseif 6==ifreqs, freqlolimit = 30; freqhilimit = 80; 
            elseif 7==ifreqs, freqlolimit = 30; freqhilimit = 50;
            elseif 8==ifreqs, freqlolimit = 50; freqhilimit = 80;              
            end
            
            % remove nirs freqs
            if 0
                for klm = 1:22,
                    f1 = klm.*Fsnirs-.4; f2 = klm.*Fsnirs+.4;
                    indxs = find( f1<Freq & Freq<=f2 );
                    Pow(indxs) = 0;
                end
            end
            
            %EEG power  ( time x chans x freqbands )
            indxs = find( freqlolimit<Freq & Freq<=freqhilimit );
            totalbandpower = sum( Pow(indxs,:), 1 );
            xxeegfclipbandpower(:,kchan,ifreqs) = totalbandpower;
            
%             %EEG power resampled to NIRS time points ( time x chans x freqbands )
%             FF = griddedInterpolant( tteegclip, xxeegfclipbandpower(:,kchan,ifreqs), 'spline','none' );
%             xxeegfclipbandpowerds(:,kchan,ifreqs) = FF(ttnirsclip);
        end
    end
    
    %% assign info data
    infocellseeg = cell(22,1);
    infocellsnirs = cell(22,1);
    infocells = cell(4,1);
    
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
    infocells{3} = tw;
    infocells{4} = overlap;
    %% assign return data
    datacells{count,1} = infocells;
    % synced data
    datacells{count,2} = tteegclip;
    datacells{count,3} = xxeegfclip;
    datacells{count,4} = ttnirsclip;
    datacells{count,5} = xxnirsfclip;
    datacells{count,6} = xxnirsfRclip;
    % eeg resampled to nirs time points
    %datacells{count,7} = xxeegfclipds;
    % eeg band power
    datacells{count,8} = xxeegfclipbandpower;
    % eeg band power resampled to nirs time points
    %datacells{count,9} = xxeegfclipbandpowerds;
    % event times
    % data clipped at first event, thus first event time is 0
    eventtimes = tevteeg - tevteeg(1);
    datacells{count,10} = eventtimes;

    count = count + 1;
    
end

