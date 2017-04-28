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
    
    experimentid = [704,804,1004,1102,1112];
    experimentid = [703,803,1003,1103,1113];
    experimentid = [703,704,803,804,1003,1004,1102,1103,1112,1113,1122,1123,1132,1133,1134];
    
    experimentid = [1132];
    
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

fignum = 380;

ifreqband = 6;
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
%% assign analysis params
flag_nirs_lockedto_eegfluctuations = 0;
flag_nirs_eeg_lagged_correlation = 0;
flag_nirs_lockedto_experimentalblocks = 1;

maxlagtime = 30;
taskduration = 30;

if 1
    % Rest Task Rest Control
    maxlagtime_blockavg = 120;
    samplenum_begin_block = [1 5 9];
elseif 0
    % Jaw Rest Smile Rest KnitBrow Rest
    maxlagtime_blockavg = 180;
    samplenum_begin_block = [1 7];
else
    error();
end
%%
figspath = '.\figs\';
figname = 'fig_eegnirscorr_';
chanstoinclude = { 'F7' 'F8'  };
chanstoinclude = { 'FZ' 'CZ' 'T4' 'F7' 'F8' };
chanstoinclude = { 'FZ' 'CZ'  };
chanstoinclude = {'FP1' 'FP2' };% 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6'}; % all but T3T4F7F8
chanstoinclude = { 'F7' 'T4'  }; % least noisy
chanstoinclude = { 'C3' 'C4' };
chanstoinclude = {'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6'}; % all but T3T4F7F8 but P
chanstoinclude = {  'P3' 'P4' 'PZ'  'O1'  'O2' }; % occipital and parietal
chanstoinclude = {'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6'}; % full 16
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6' 'P3' 'P4' 'PZ' }; % full 19
chanstoinclude = { 'T3SUPER' 'T4SUPER' 'F7SUPER' 'F8SUPER' }; %
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'T3SUPER' 'T4SUPER' 'F7SUPER' 'F8SUPER'}; %
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' }; %
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6' 'P3' 'P4' 'PZ'  'T3SUPER' 'T4SUPER' 'F7SUPER' 'F8SUPER'}; % full 19


%%
if flag_nirs_lockedto_experimentalblocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  EEG BAND POWER - NIRS
    %  NIRS LOCKED TO EXPERIMENTAL BLOCKS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    %% loop over files
    nirsOeventavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
    nirsReventavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
    nirsOblockavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
    nirsRblockavgbar = [];%zeros(nlags.*2+1, 1); % each experiment gets filled asa col
    for experimentnun = [ 1 ],
        datacells = datacellsarray{experimentnun};
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
        if 1,
            length(eegchanlabels)
            eegchanlabels
            length(nirschanlabels)
            nirschanlabels
            ttnirsclip(end)./60            
            return;
        end
        %% EMG
        % freq band
        if 0
            
            figure(23211); clf; count = 1;
            for jj = 1:numchaneeg, % EMG
                chan2 = upper ( char( eegchanlabels(jj) ) );
                if [strfind( 'LOC', chan2 ) ...
                        strfind( 'CH24', chan2 ) ...
                        strfind( 'ROC', chan2 ) ...
                        strfind( 'CH25', chan2 )] ...
                        fprintf('[%s]\n',chan2 );
                    eegdata = ( xxeeg(:,jj,ifreqband) );
                    inddxx=[1:floor(Fsnirs.*12.*60)]; hold on; grid off;
                    inddxx=[1:length(eegdata)];
                    subplot(4,1,count),
                    count = count + 1;
                    % event markers
                    hold on;
                    gry2 = [.7 .9 .7]; ymin = min(eegdata); ymax = max(eegdata); eventlinewidth = 2;for jjkk=1:length(eventtimes), plot( [eventtimes(jjkk) eventtimes(jjkk)], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth); end
                    %plot( ttnirsclip(inddxx),xxeeg(inddxx,jj,ifreqband));
                    plot(ttnirsclip(inddxx),eegdata(inddxx),'k')
                    title([ chan2 ' ' fpath ' ' fpath2 ],'interpreter','none');
                    ymin = min(eegdata); ymax = max(eegdata);
                    % ymin = min(eegdata); ymax = 20.*mean(eegdata);
                    %ymin = -1;  ymax = 1000;
                    %xmin = 0; xmax = 870;
                    axis([ttnirsclip(1) ttnirsclip(end) ymin ymax])
                    pause(.1);
                    % retval = input('Press <Enter> to continue...');
                end
            end
            
            if 0,
                fnamefig =  ['fig_EMG_' num2str(experimentid) '.tiff'];
                set( gcf, 'PaperPositionMode', 'auto')
                print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
            end
        end
        %% EMG
        % spectrogram
        if 1
            figure(23213); clf; count = 1;
            for jj = 1:numchaneeg, % EMG
                chan2 = upper ( char( eegchanlabels(jj) ) );
                if [strfind( 'LOC', chan2 ) ...
                        strfind( 'CH24', chan2 ) ...
                        strfind( 'ROC', chan2 ) ...
                        strfind( 'CH25', chan2 )] ...
                        fprintf('[%s]\n',chan2 );
                    eegdata = ( xxeegorig(:,jj) );
                    windowspec = ceil(Fs.*1); nooverlapspec = ceil(windowspec./2); 
                    deltafreq = Fs./(2.*windowspec); F = [0:deltafreq:80];
                    [S,F,T,P] = spectrogram(eegdata, windowspec, nooverlapspec, F, Fs);
                    % S: Fourier coeffs, F: freqs, T: times, P: PSD
                    subplot(4,1,count),
                    count = count + 1;
                    surf(T,F,10.*log10(P),'edgecolor','none'); view(0,90); %xlabel('Time (s)'); 
                    ylabel('Hz'); title(['EMG spectrogram ' chan2 ]);
                    set(gca,'Xtick',[0:taskduration:T(end)]); set(gca,'TickDir','out');
                    axis([T(1) T(end) 0 80]);
                    drawnow;
                    pause(.1);
                    % retval = input('Press <Enter> to continue...');
                end
            end
            if 0,
                fnamefig =  ['fig_EMG_' num2str(experimentid) '.tiff'];
                set( gcf, 'PaperPositionMode', 'auto')
                print(  '-dtiff', '-zbuffer', '-r600',  [figspath  fnamefig] );
            end
        end
        %% EEGPOWER-NIRS
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
                chan2 = upper ( char( eegchanlabels(jj) ) );
                
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
                        %  ymax = 1000;
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
                     %   ymin = min(eegdata); ymax = max(eegdata);
                     %   ymax = 100;
                        ymin = min(eegdata); ymax = 8.*mean(eegdata);                        
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
                    
                    %% AVERAGED NIRS SIGNAL
                    % NIRS LOCKED TO EXPERIMENT BLOCK
                    bucket1 = []; bucket2 = [];
                    for inn = samplenum_begin_block,
                        indxx = [ floor(eventtimes(inn).*Fsnirs) + 1 : floor(eventtimes(inn).*Fsnirs) + nlags_blockavg];
                        if 1<=indxx(1)&indxx(end)<=length(ttnirsclip),
                            bucket1 = [bucket1 nirsOdata(indxx)];
                            bucket2 = [bucket2 nirsRdata(indxx)];
                        end
                    end
                    nirsOblockavg = mean(bucket1,2); % mean over events for one experiment one chan
                    nirsRblockavg = mean(bucket2,2); % mean over events for one experiment one chan
                    
                    if 1, % plot experiment block avg
                        figure(2024); clf; hold on;
                        title( { [ chan1]; [fpath ' ' fpath2]; ['AVG EXPERIM. BLOCK']},'interpreter','none','fontsize',12);
                        ymax = max([nirsOblockavg; nirsRblockavg]); ymin = min([nirsOblockavg; nirsRblockavg]);
                                 
                        mytim_blockavg = [1:length(nirsOblockavg)]./Fsnirs;                        

                        plot([mytim_blockavg(1) mytim_blockavg(end)],[0 0],'linestyle','-','color',gry2,'linewidth', eventlinewidth);
                        for iijjkk = 1:5, deltatime = iijjkk.*30; plot( [deltatime deltatime], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth);                        end
                      
                        ns = ceil( 10.*Fsnirs );
                        plot(  mytim_blockavg, ( mysmooth( nirsOblockavg, ns) ), 'r' );
                        plot(  mytim_blockavg, ( mysmooth( nirsRblockavg, ns) ), 'b' );
                        %axis([0 mytim_blockavg(end) ymin ymax]);
                        
                        pause(.1);
                        %  retval = input('Press <Enter> to continue...');
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
                    
                    %% RESAMPLE TO GLOBAL SAMPLE RATE (INDEPT OF EXPERIM)
                    % NIRS LOCKED TO EXPERIMENT BLOCK  
                    mytim_blockavg = [1:length(nirsOblockavg)]./Fsnirs;
                    nlags0_blockavg = ceil( Fsnirs0.*maxlagtime_blockavg );
                    mytim0_blockavg = [1:(nlags0_blockavg+1)]./Fsnirs0;
                    
                    FF = griddedInterpolant( mytim_blockavg, nirsOblockavg, 'spline','none' );
                    nirsOblockavg0 = FF(mytim0_blockavg)';
                    FF = griddedInterpolant( mytim_blockavg, nirsRblockavg, 'spline','none' );
                    nirsRblockavg0 = FF(mytim0_blockavg)';
                    
                    if 1, % plot experiment block avg
                        figure(1024); clf; hold on;
                        title( { [ chan1]; [fpath ' ' fpath2]; ['AVG EXPERIM. BLOCK']},'interpreter','none','fontsize',12);
                        ymax = max([nirsOblockavg0; nirsRblockavg0]); ymin = min([nirsOblockavg0; nirsRblockavg0]);
                        
                        % mark the experimental blocks
                        gry2 = [.7 .9 .7]; eventlinewidth = 2;
                        plot([mytim0_blockavg(1) mytim0_blockavg(end)],[0 0],'linestyle','-','color',gry2,'linewidth', eventlinewidth);
                        for iijjkk = 1:5, deltatime = iijjkk.*30; plot( [deltatime deltatime], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth);                        end
                        
                        ns = 1;%ceil( 1.*Fsnirs0);
                        plot( mytim0_blockavg, ( mysmooth( nirsOblockavg0, ns) ), 'r' );
                        plot( mytim0_blockavg, ( mysmooth( nirsRblockavg0, ns) ), 'b' );
                        %axis([0 mytim0_blockavg(end) ymin ymax]);
                        
                        % axis([0 maxlagtime_blockavg ymin ymax]);
                        pause(.1);
                        %  retval = input('Press <Enter> to continue...');
                    end
                    nirsOblockavgbar = [nirsOblockavgbar nirsOblockavg0]; % this accumulates each chan as a column
                    nirsRblockavgbar = [nirsRblockavgbar nirsRblockavg0]; % this accumulates each chan as a column

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
    
    %% inter experiment average
    % NIRS LOCKED TO EXPERIMENT BLOCK
    nirsOblockavgbar0 = mean( nirsOblockavgbar, 2 );
    nirsRblockavgbar0 = mean( nirsRblockavgbar, 2 );
    if 1
        
        ns = ceil( 10.*Fsnirs0);
        
        % HBO
        figure(29221); clf; hold on; %title([ fpath ' ' fpath2 ],'interpreter','none'); grid on;
        %ymax = max(max([nirsOblockavgbar])); ymin = min(min([nirsOblockavgbar]));
        
        % mark the experimental blocks
        gry2 = [.7 .9 .7]; eventlinewidth = 2; plot([mytim0_blockavg(1) mytim0_blockavg(end)],[0 0],'linestyle','-','color',gry2,'linewidth', eventlinewidth);
        for iijjkk = 1:5, deltatime = iijjkk.*30; plot( [deltatime deltatime], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth);                        end
        
        plot( mytim0_blockavg, nirsOblockavgbar, 'r','linewidth',1 );
        plot( mytim0_blockavg, mysmooth(nirsOblockavgbar0,ns), 'r','linewidth',4 );
        ymin = -5e-5; ymax = 5e-5;
        axis([0 maxlagtime_blockavg ymin ymax]);
        
        % HBR
        figure(29222); clf; hold on; title([ fpath ' ' fpath2 ],'interpreter','none');
        %ymax = max(max([nirsRblockavgbar])); ymin = min(min([nirsRblockavgbar]));
        
        % mark the experimental blocks
        gry2 = [.7 .9 .7]; eventlinewidth = 2; plot([mytim0_blockavg(1) mytim0_blockavg(end)],[0 0],'linestyle','-','color',gry2,'linewidth', eventlinewidth);
        for iijjkk = 1:5, deltatime = iijjkk.*30; plot( [deltatime deltatime], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth);                        end
        
        plot( mytim0_blockavg, nirsRblockavgbar, 'b','linewidth',1 );
        plot( mytim0_blockavg, mysmooth(nirsRblockavgbar0,ns), 'b','linewidth',4 );
        ymin = -1e-5; ymax = 1e-5;
        axis([0 maxlagtime_blockavg ymin ymax]);
        
        % HBO HBR
        figure(29224); clf; hold on; title([ fpath ' ' fpath2 ],'interpreter','none');
        %ymax = max(max([nirsOblockavgbar0; nirsRblockavgbar0])); ymin = min(min([nirsOblockavgbar0; nirsRblockavgbar0]));
        
        % mark the experimental blocks
        gry2 = [.7 .9 .7]; eventlinewidth = 2; plot([mytim0_blockavg(1) mytim0_blockavg(end)],[0 0],'linestyle','-','color',gry2,'linewidth', eventlinewidth);
        for iijjkk = 1:5, deltatime = iijjkk.*30; plot( [deltatime deltatime], [ymin ymax],'linestyle','-','color',gry2,'linewidth', eventlinewidth);                        end
        
        plot( mytim0_blockavg, nirsOblockavgbar0, 'r','linewidth',4 );
        plot( mytim0_blockavg, nirsRblockavgbar0, 'b','linewidth',4 );
        ymin = -5e-5; ymax = 5e-5;
        axis([0 maxlagtime_blockavg ymin ymax]);
        pause(.1);
    end
    return;
end


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
        for experimentnun = [ 1 ],
            datacells = datacellsarray{experimentnun};
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
        for experimentnun = [ 1 ],
            datacells = datacellsarray{experimentnun};
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
