function datacells = geteegnirsdata_mixedsourcedetectors( ...
    experimentid,...
    dtbegeeg,...
    dtendeeg,...
    filterlofreqeeg,...
    filterhifreqeeg,...
    filternotcheeg,...
    numsampleswindoweegpower,...
    filterlofreqnirs,...
    filterhifreqnirs ,...
    filternotchnirs...
    )

 
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

