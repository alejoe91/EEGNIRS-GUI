clear; format compact; fprintf('\n');
%% DATA FILE PARAMS

% output path
outputpath = 'C:\Users\aomurtag\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\corrmat\';

% channel names and locations
fchans =  'positions_wholehead_flat.txt';

% Declare empty ARRAYS which hold the files and their properties
fnames = {};
fpaths = {};
badchanscellarr = {};

% Fill in the arrays for input files
if 0
    indxfile = 1;
    fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
    fname = ['Omurtag327eyescloseonlyeeg_20140327_161725'];
    fnames = [fnames; fname]; fpaths = [fpaths; fpath];
    badchans = {}; % example: badchans = {'FP1', 'T6' 'F4' 'FP1' 'F1' 'F8' 'P3' 'EKG1' };
    badchanscellarr{indxfile} = badchans;
end

if 0
    indxfile = 1;
    fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
    fname = ['Omurtag327restingstateboth_20140327_162652'];
    fnames = [fnames; fname]; fpaths = [fpaths; fpath];
    badchans = {};
    badchanscellarr{indxfile} = badchans;
    
    indxfile = 2;
    fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
    fname = ['Omurtag327restingstateeyesopenboth_20140327_163413'];
    fnames = [fnames; fname]; fpaths = [fpaths; fpath];
    badchans = {};
    badchanscellarr{indxfile} = badchans;
end

if 0
    indxfile = 1;
    fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
    fname = ['Omurtag327taskeyesopenboth_20140327_165145'];
    fnames = [fnames; fname]; fpaths = [fpaths; fpath];
    badchans = {};
    badchanscellarr{indxfile} = badchans;
    
    indxfile = 2;
    fpath = 'F:\LabData\LabData\AhmetOmurtag0327\EEG\';
    fname = ['Omurtag327taskeyescloseboth_20140327_164238'];
    fnames = [fnames; fname]; fpaths = [fpaths; fpath];
    badchans = {};
    badchanscellarr{indxfile} = badchans;
end

if 1  
    indxfile = 1;
    fpath = 'F:\LabData\LabData\AhmetOmurtag0416\EEG\';
    fname = ['AhmetOmurtag41615minutestask_20140416_191933'];
    fnames = [fnames; fname]; fpaths = [fpaths; fpath];
    badchans = {};
    badchanscellarr{indxfile} = badchans;  
end


%% ANALYSIS params
Tw = 1; % analysis window (sec)
Nwmax = 100000000; % max num windows to use (large number implies use max available)
nlags = 40; % max number of lags to compute (samples)

%% PREPROCESS PARAMS
% data clip
dtbeg = 60; % sec
dtend = 60; % sec
% Bandpass filter and notch (normal bandpass range: 0.5-70 Hz)
filterlofreq = .5; % .5; % Hz
filterhifreq = 70; % 70; % Hz
filternotch = 60; % 0: no notch;
flag_normalize_datawindow = true;

% plot lagged correlation params
flag_plotlaggedcorr = 1;

%% read and analyze all files in a loop
for ii = 1:length(fnames),
    
    fname = char( fnames{ii} );
    fpath = char( fpaths{ii} );
    fullfname = [ fpath  fname '.edf' ];
    badchans = badchanscellarr{ii};
    fprintf('[%s]\n', fullfname);
    badchans
    
    %% READ DATA
    numrecordstoread = 10000; % to read all records set numrecordstoread=[]
    numchanmax = 26; %21;
    [genheadr1, sigheadr1, totnumrecs1, xx, Fs, numchan] = readEDFfile2( fullfname, numrecordstoread, numchanmax );
    
    %% CLIP BEGINNING AND END
    nbeg = Fs.*dtbeg;
    nend = Fs.*dtend;
    xx = xx( nbeg:end-nend, : );
    %% FILTER
    xxf = mybutterandiirnotchfilters( xx, [filterlofreq filterhifreq filternotch 1], 6, Fs );
    %xxf = mybandstopfilter( xxf, 6, 7, 6, Fs );
    %% initialize the names and positions of predetermined full set of 128 channels
    [ labels, xchan ] = inputchanlocations( fchans );
    subsetindxs = [];
    fprintf('Channels Used: \n');
    for ii = 1:numchan,
        chan1 = sigheadr1(ii).label;
        [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
        if ifound1,
            subsetindxs = [subsetindxs; ii];
            fprintf( '%d %s, ', ii, chan1);
        end
    end
    fprintf('\n');
    %% VIEW SIG
    if 0
        nbeg = 1; %Fs.*1;
        nend = nbeg + Fs.*10;
        tt = [nbeg:nend]./Fs;
        for ii = 1:numchan,
            chan1 = sigheadr1(ii).label; [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
            if ifound1,
                %            figure(110); clf; plot( tt, xxf(nbeg:nend, ii) ); title( [ num2str(ii) ':   ' chan1 ] );
                figure(1102); clf; pmtm( xxf(nbeg:nend, ii), 3.5, [], Fs); title( [ num2str(ii) ':   ' chan1 ] );
                pause(.5);
            else
                figure(112); clf; plot( tt, xxf(nbeg:nend, ii) ); title( [ num2str(ii) ':   ' chan1 ] );
                pause(.5);
            end
        end
    end
    %% CLEAN UP (eg EYEBLINKS)
    if 0
        xtmp = xxfnc( : , subsetindxs);
        xoc1 = xxfnc( : , 23 );
        xoc2 = xxfnc( : , 25 );
        [xtmpc] = remove_ocular_adapt_filter( xtmp, xoc1, xoc2, Fs ); % Required DSP toolbox
        xxfnc( : , subsetindxs ) = xtmpc;
    else
        xxfnc = xxf;
    end
    %% MAX OF THE LAGGED CORR WITHIN WINDOW
    numsampW = Fs.*Tw;
    nstds = [];
    
    maxdim = min( Nwmax, floor( length(xxfnc)./numsampW ) );
    
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
    while iw<=Nwmax & iend<=length(xxfnc),
        
        if 0==mod(iw,10), fprintf('%d ',iw); end;
        
        dat = xxfnc( ibeg:iend, : );
        if flag_normalize_datawindow, dat = mynormalize( dat ); end
        
        for ii = 1:numchan,
            chan1 = sigheadr1(ii).label;
            %            [ xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
            if isempty( intersect(subsetindxs, ii) ), ifound1=0; else, ifound1=1; end;
            indxbadchan1 = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
            
            for jj = ii:numchan,
                chan2 = sigheadr1(jj).label;
                %                [ xxx2, iii2 , ifound2 ] = getchanlocation( chan2, labels, xchan );
                if isempty( intersect(subsetindxs, jj) ), ifound2=0; else, ifound2=1; end;
                indxbadchan2 = findstrincellarray( badchans, upper(chan2) );
                
                %                fprintf('[%s] [%s]\n',chan1,chan2);
                
                if ifound1 & ifound2 & 0==indxbadchan1 & 0==indxbadchan2,
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
                    
                    if flag_plotlaggedcorr, % plot correlation vs. lag
                        figure(880); clf;
                        plot( Lags./Fs, XCF, 'k' ); hold on;
                        maxlim = max( Lags./Fs); axis([-maxlim maxlim -1 1]); grid on;
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
        iw = iw + 1;
        % bounds of data for next window
        ibeg = (iw-1).*numsampW + 1;
        iend = ibeg + numsampW - 1;
    end
    
    Nw = iw - 1; % Total num of windows
    fprintf('\n');
    %% SAVE corr data
    savestring = ['save ' outputpath 'corrmat_freqs_' num2str(filterlofreq) '-' num2str(filterhifreq) '-Tw'  num2str(Tw) '_' fname '.mat sigheadr1 Tw Nw Fs badchans labels xchan subsetindxs corrmax lagcorrmax'];
    fprintf('[%s]\n', savestring);
    eval( savestring );
    
end % end read and analyze each files in a loop

