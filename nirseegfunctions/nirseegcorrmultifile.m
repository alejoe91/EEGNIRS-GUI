%clear;
%format compact;
fprintf('\n');
%% compute get all data
%datasets = [703,704,803,804,903,904,1003,1004];
%datasets = [1003,1004];
%datacells = geteegnirsdata( datasets );
%% load data
matlabdatafolder = 'D:\matlab\matlabdata';
%{
loadstring = ['load ' matlabdatafolder '\nirseegcorrmultfile_datacells1.mat;' ];
fprintf('[%s]\n', loadstring);
eval( loadstring );
%}
%{
    % save data
    matlabdatafolder = 'D:\matlab\matlabdata';
    savestring = ['save '  matlabdatafolder '\nirseegcorrmultfile_datacells1.mat ' ' datacells'];
    fprintf('[%s]\n', savestring);
    eval( savestring );
%}
%{
    %% assign info data
    infocellseeg = cell(8,1);
    infocellsnirs = cell(8,1);
    infocells = cell(3,1);
    infocellsnirs{1} = fpath;
    infocellsnirs{2} = fpath2;
    infocellsnirs{3} = fname;
    infocellseeg{4} = fpatheeg;
    infocellseeg{5} = fnameeeg;
    infocellsnirs{6} = badchans;
    infocellseeg{7} = Fs;
    infocellsnirs{7} = Fsnirs;
    infocells{1} = infocellseeg;
    infocells{2} = infocellsnirs;
    %% assign  data
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
%}
%%
if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  EEG - NIRS lAGGED CORR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% assign analysis params
    idisplaytitle = 1;
    maxlagtime = 40;
    figspath = '.\figs\';
    figname = 'fig_eegnirscorr_';
    flag_useprefrontalonly = 1;
    %% loop over files
    xcorrelbar = [];%zeros(nlags.*2+1, 1);
    xcorrelRbar = [];%zeros(nlags.*2+1, 1);
    xcorrel3bar = [];%zeros(nlags.*2+1, 1);
    for filenum = [1 3 5],
        
        %% retrieve data
        infocells = datacells{filenum,1}; infocellseeg = infocells{1}; Fs = infocellseeg{7}; eegchanlabels = infocellseeg{9}; numchaneeg = infocellseeg{8};;
        xxeegfclipds = datacells{filenum,7}; % EEG resampled
        
        infocellsnirs = infocells{2}; Fsnirs = infocellsnirs{7}; nirschanlabels = infocellsnirs{9}; numchannirs = infocellsnirs{8}; xxnirsfclip = datacells{filenum,5};
        xxnirsfRclip = datacells{filenum,6};
        badchans = infocellsnirs(6);
        %%
        fpath = char( infocellsnirs{1} ); fpath2 = char( infocellsnirs{2} );
        fprintf('%d [%s]\n',filenum,[fpath fpath2]);
        %% Cconvert from MCN system to international 10-20:  T3, T4, T5 and T6— instead of T7, T8, P7 and P8
        eegchanlabels = convertMCNto1020(eegchanlabels);
        nirschanlabels = convertMCNto1020(nirschanlabels);
        %% Lagged corr EEG-NIRS
        nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
        
        for ii = 1:numchannirs, % NIRS
            chan1 = upper( char( nirschanlabels(ii) ) );
            
            indxbadchan = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
            if 0~=indxbadchan, continue; end
            
            if flag_useprefrontalonly&~isprefrontal( char(chan1) ), continue; end
            
            [ xcorrel3, Lags ] = xcorr(   xxnirsfRclip(:,ii), xxnirsfclip(:,ii), nlags, 'coeff' );
            xcorrel3bar = [xcorrel3bar xcorrel3];
            
            %figure(1020);clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel3,'k');
            
            for jj = 1:numchaneeg, % EEG
                chan2 = upper ( char( eegchanlabels(jj) ) );
                
                indxbadchan = findstrincellarray( badchans, upper(chan2) ); % return zero if chan1 is NOT a bad channel
                if 0~=indxbadchan, continue; end
                
                if strcmp( chan1, chan2 ),
                    % fprintf('[%s][%s]\n',chan1,chan2 );
                    
                    [ xcorrel, Lags ] = xcorr(  xxnirsfclip(:,ii), xxeegfclipds(:,jj),  nlags, 'coeff' );
                    [ xcorrelR, Lags ] = xcorr(   xxnirsfRclip(:,ii), xxeegfclipds(:,jj), nlags, 'coeff' );
                    xcorrelbar = [xcorrelbar xcorrel];
                    xcorrelRbar = [xcorrelRbar xcorrelR];
                    
                    %figure(1020); hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel,'r'); plot(Lags./Fsnirs,xcorrelR,'b');
                    % figure(1010); clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel,'r'); plot(Lags./Fsnirs,xcorrelR,'b');
                    %pause(.1);
                end
            end
        end
        
    end
    xcorrelbar0 = mean(xcorrelbar,2);
    xcorrelRbar0 = mean(xcorrelRbar,2);
    xcorrel3bar0 = mean(xcorrel3bar,2);
    figure(2040);clf; hold on; plot( xcorrelbar ); plot( xcorrelbar0,'linewidth',11 )
    figure(2010);clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrelbar0,'r'); plot(Lags./Fsnirs,xcorrelRbar0,'b');
    figure(2020);clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel3bar0,'k');
end
%
if 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  EEG BAND POWER - NIRS lAGGED CORR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% assign analysis params
    idisplaytitle = 1;
    maxlagtime = 50;
    figspath = '.\figs\';
    figname = 'fig_eegnirscorr_';
    flag_useprefrontalonly = 0;
    chanstoinclude = { 'F7' 'F8'  };
    chanstoinclude = { 'FZ' 'CZ' 'T4' 'F7' 'F8' };
    chanstoinclude = { 'FZ' 'CZ'  };
    chanstoinclude = { 'T3' 'T4' 'F7' 'F8' };
    chanstoinclude = { 'O1' 'O2' };
    %% loop over files
    xcorrelbar = [];%zeros(nlags.*2+1, 1);
    xcorrelRbar = [];%zeros(nlags.*2+1, 1);
    xcorrel3bar = [];%zeros(nlags.*2+1, 1);
    for filenum = [2 4 6 8],
        
        %% retrieve data
        infocells = datacells{filenum,1}; infocellseeg = infocells{1}; Fs = infocellseeg{7}; eegchanlabels = infocellseeg{9}; numchaneeg = infocellseeg{8};;
        xxeegfclipbandpowerds = datacells{filenum,9}; % EEG bandpower resamples
        
        infocellsnirs = infocells{2}; Fsnirs = infocellsnirs{7}; nirschanlabels = infocellsnirs{9}; numchannirs = infocellsnirs{8};
        xxnirsfclip = datacells{filenum,5}; xxnirsfRclip = datacells{filenum,6};
        badchans = infocellsnirs(6);
        %%
        fpath = char( infocellsnirs{1} ); fpath2 = char( infocellsnirs{2} );
        fprintf('%d [%s]\n',filenum,[fpath fpath2]);
        %% Cconvert from MCN system to international 10-20:  T3, T4, T5 and T6— instead of T7, T8, P7 and P8
        eegchanlabels = convertMCNto1020(eegchanlabels);
        nirschanlabels = convertMCNto1020(nirschanlabels);
        %% Lagged corr EEGPOWER-NIRS
        nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
        
        ifreqband = 6; %6
        
        for ii = 1:numchannirs, % NIRS
            chan1 = upper( char( nirschanlabels(ii) ) );
            
            indxbadchan = findstrincellarray( badchans, upper(chan1) ); % return zero if chan1 is NOT a bad channel
            if 0~=indxbadchan, continue; end
            
            indxchanstoinclude = findstrincellarray( chanstoinclude, upper(chan1) ); % return zero if chan1 is NOT a bad channel
            if 0==indxchanstoinclude, continue; end
            
            if flag_useprefrontalonly&~isprefrontal( char(chan1) ), continue; end
            
            [ xcorrel3, Lags ] = xcorr(   xxnirsfRclip(:,ii), xxnirsfclip(:,ii), nlags, 'coeff' );
            xcorrel3bar = [xcorrel3bar xcorrel3];
            
            %figure(1020);clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel3,'k');
            
            for jj = 1:numchaneeg, % EEG
                chan2 = upper ( char( eegchanlabels(jj) ) );
                
                if strcmp( chan1, chan2 ),
                     fprintf('[%s][%s]\n',chan1,chan2 );
                    
                    [ xcorrel, Lags ] = xcorr(  xxnirsfclip(:,ii), xxeegfclipbandpowerds(:,jj,ifreqband),  nlags, 'coeff' );
                    [ xcorrelR, Lags ] = xcorr(   xxnirsfRclip(:,ii), xxeegfclipbandpowerds(:,jj,ifreqband), nlags, 'coeff' );
                    xcorrelbar = [xcorrelbar xcorrel];
                    xcorrelRbar = [xcorrelRbar xcorrelR];
                    
                    figure(1020); clf; hold on; title([ chan1 ' ' fpath ' ' fpath2 ],'interpreter','none'); grid on; 
                    plot(Lags./Fsnirs,xcorrel,'r'); plot(Lags./Fsnirs,xcorrelR,'b'); 
                    axis([-maxlagtime maxlagtime -.2 .2])
                    % figure(1010); clf; hold on; title([ chan1 ]); grid on; plot(Lags./Fsnirs,xcorrel,'r'); plot(Lags./Fsnirs,xcorrelR,'b');
                   pause(.1);
                end
            end
        end
        
    end
    xcorrelbar0 = mean(xcorrelbar,2);
    xcorrelRbar0 = mean(xcorrelRbar,2);
    xcorrel3bar0 = mean(xcorrel3bar,2);
    figure(3040);clf; hold on; plot( Lags./Fsnirs, xcorrelbar ,'r' ); plot( Lags./Fsnirs, xcorrelbar0,'linewidth',5,'color','r' )
    figure(3050);clf; hold on; plot(Lags./Fsnirs, xcorrelRbar, 'b' ); plot( Lags./Fsnirs,xcorrelRbar0,'linewidth',5,'color','b' )
    %%
    gryr = [.99 .5 .5];    gryb = [.5 .5 .99];

        figure(3060);clf; hold on; offst = -.3;
        plot( Lags./Fsnirs, xcorrelbar ,'color',gryr,'linewidth',.1,'linestyle','-' ); plot( Lags./Fsnirs, xcorrelbar0,'linewidth',5,'color','r' )
plot(Lags./Fsnirs, xcorrelRbar + offst ,'color',gryb,'linewidth',.1,'linestyle','-' ); plot( Lags./Fsnirs,xcorrelRbar0+ offst,'linewidth',5,'color','b' )
    %%
    
    figure(3010);clf; hold on; grid on; plot(Lags./Fsnirs,xcorrelbar0,'r'); plot(Lags./Fsnirs,xcorrelRbar0,'b');
    figure(3020);clf; hold on; grid on; plot(Lags./Fsnirs,xcorrel3bar0,'k');
end