format compact;
fprintf('\n');
flag_readdatafromfile = 1;
%% get data
%% compute and store
if flag_readdatafromfile,
    clear;
    dtbegeeg = 10; %sec
    dtendeeg = 10;
    filterlofreqeeg = .5;
    filterhifreqeeg = 70;
    filternotcheeg = 60;
    numsampleswindoweegpower = 50;
    filterlofreqnirs = 0.01;
    filterhifreqnirs = 0.08;
    filternotchnirs = 0;
      infostring = [ '_eegfilt' num2str(filterlofreqeeg) '-' num2str(filterhifreqeeg) '_nirsfilt' num2str(filterlofreqnirs) '-' num2str(filterhifreqnirs) ];
   % infostring = [ '_eegfilt' num2str(filterlofreqeeg) '-' num2str(filterhifreqeeg) '_nirsfilt' num2str(filterlofreqnirs) '-' num2str(filterhifreqnirs) '_eegpowwindow' num2str(numsampleswindoweegpower) ];
    experimentid = [703,704,803,804,1003,1004,1102];
    experimentid = [1004];
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

ifreqband = 6; %6: 30-70
%{            
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
%}
%% assign analysis params
maxlagtime = 100;
figspath = '.\figs\';
figname = 'fig_eegnirs_mapcorr_';

chanstoinclude = { 'F7' 'F8'  };
chanstoinclude = { 'FZ' 'CZ' 'T4' 'F7' 'F8' };
chanstoinclude = { 'FZ' 'CZ'  };
chanstoinclue = { 'FZ' 'CZ' };
chanstoinclude = {'FP1' 'FP2' };% 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6'}; % all but T3T4F7F8
chanstoinclude = { 'F7' 'T4'  }; % least noisy
chanstoinclude = {'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6'}; % all but T3T4F7F8
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' }; % solid
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6' 'P3' 'P4' 'PZ' }; % full 19
chanstoinclude = { 'T3' 'T4' 'F7' 'F8' 'FP1' 'F3' 'T5' 'C3' 'FZ' 'O1' 'FP2'  'CZ' 'F4' 'C4' 'O2' 'T6'}; % full 16
%% initialize the names and positions of predetermined full set of 128 channels
[ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
labels = convertMCNto1020(labels);
labels = charcellarraytoupper( labels );
covsize = length(chanstoinclude);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EEG BAND POWER PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kminutes = 12:12,
    % define covariance matrices
    cove = zeros(covsize); covno = zeros(covsize); covnr = zeros(covsize);
    % define concatenated time series
    xxe = []; xxno = []; xxnr = []; numsamples = [];
    %% loop over files
    for experimentnum = [ 1 ],
        datacells = datacellsarray{experimentnum};
        filenum = 1;
        %% retrieve data
        infocells = datacells{filenum,1}; infocellseeg = infocells{1}; Fs = infocellseeg{7}; eegchanlabels = infocellseeg{9}; numchaneeg = infocellseeg{8};;
        xxeeg0 = datacells{filenum,9}; % EEG bandpower resampled
        xxeeg = xxeeg0(:,:,ifreqband); clear xxeeg0;
        infocellsnirs = infocells{2}; Fsnirs = infocellsnirs{7}; nirschanlabels = infocellsnirs{9}; numchannirs = infocellsnirs{8};
        xxnirsfclip = datacells{filenum,5}; xxnirsfRclip = datacells{filenum,6};
        %% char cell array to upper
        eegchanlabels = charcellarraytoupper( eegchanlabels );
        nirschanlabels = charcellarraytoupper( nirschanlabels );
        %% take  part of data - to study convergence
        numsamptouse = floor(kminutes.*60.*Fsnirs); % floor(numsamptot./1);
        xxeeg = xxeeg(1:numsamptouse,:,:);
        xxnirsfclip = xxnirsfclip(1:numsamptouse,:);
        xxnirsfRclip = xxnirsfRclip(1:numsamptouse,:);
        % num samples
        numsamples = [numsamples; length(xxeeg)];
        %%
        badchans = infocellsnirs(6);
        %% screen display progress
        fpath = char( infocellsnirs{1} ); fpath2 = char( infocellsnirs{2} );
        fprintf('%d [%s]\n',filenum,[fpath fpath2]);
        %% Cconvert from MCN system to international 10-20:  T3, T4, T5 and T6— instead of T7, T8, P7 and P8
        eegchanlabels = convertMCNto1020(eegchanlabels);
        nirschanlabels = convertMCNto1020(nirschanlabels);
        %% select subset of chans. after this step cols arranged in order of chanstoinclude
        xxeeg = selectsubsetofchans( xxeeg, eegchanlabels, chanstoinclude );
        xxnirsfclip = selectsubsetofchans( xxnirsfclip, nirschanlabels, chanstoinclude );
        xxnirsfRclip = selectsubsetofchans( xxnirsfRclip, nirschanlabels, chanstoinclude );
        %% normalize
        % xxeeg = mynormalize(xxeeg); xxnirsfclip = mynormalize(xxnirsfclip); xxnirsfRclip = mynormalize(xxnirsfRclip);
        %% normalize
        % xxeeg = meansubtract(xxeeg); xxnirsfclip = meansubtract(xxnirsfclip); xxnirsfRclip = meansubtract(xxnirsfRclip);
        %% accumulate covariance matrices
%        cove = cove + cov(xxeeg); covno = covno + cov(xxnirsfclip); covnr = covnr + cov(xxnirsfRclip);
        cove = cove + (xxeeg'*xxeeg); covno = covno + (xxnirsfclip'*xxnirsfclip); covnr = covnr + (xxnirsfRclip'*xxnirsfRclip);
        %% accumulate time series
        xxe = [xxe; xxeeg]; xxno = [xxno; xxnirsfclip]; xxnr = [xxnr; xxnirsfRclip];
    end
    %% EEG PCA
    %    [evece,aae,lame] = princomp(xxe); %    [EE,EV] = eig(cove);    lame = diag(EV);     [lame,indx] = sort(lame,1,'descend');    evece = EE(:,indx);
    [EE,EV] = eig(cove);    lame = diag(EV);     [lame,indx] = sort(lame,1,'descend');    evece = EE(:,indx);
%    [evece,lame] = pcacov(cove);
    aae = xxe * evece;
    lamenorm = lame./sum(lame);
    %% NIRS PCA
    [evecno,lamno] = pcacov(covno);
    aano = xxno * evecno;
    lamnonorm = lamno./sum(lamno);
    [evecnr,lamnr] = pcacov(covnr);
    aanr = xxnr * evecnr;
    lamnrnorm = lamnr./sum(lamnr);
    %% PLOT PCA variances
    fh = figure(10101);clf;  hold on; set(fh,'color','w');
    ppe = plot( cumsum(lamenorm),'k.-');
    ppno = plot( cumsum(lamnonorm),'r.-');
    ppnr = plot( cumsum(lamnrnorm),'b.-');
    legend([ppe ppno ppnr],'EEG gamma','HbO','HbR','location','southeast'); legend boxoff
    title('Convergence of the PC variance','fontsize',16);
    xlabel('PC#','fontsize',16); ylabel('Variance as Fraction of the Total','fontsize',16);
    set(gca,'fontsize',16);
    %plot( cumsum(lamenorm),'r.-');    plot( cumsum(lamnonorm),'r.-');
    %% PLOT MAP
    labelclr = 'k';
    if 1
        %plot6maps( 110, 'EEG', evece, lamenorm, chanstoinclude, labels, xchan );
        plot6maps_horiz( 110, 'EEG', evece, lamenorm, chanstoinclude, labels, xchan, labelclr );
        drawnow
        if 0
            figname = ['fig_eeg_freq' num2str(ifreqband) '_pcmaps'];    set( gcf, 'PaperPositionMode', 'auto');
            print(  '-dtiff', '-zbuffer', '-r600',  [figspath  figname  '.tiff'] );
        end
    end
    if 1
        plot6maps_horiz( 120, 'HbO', evecno, lamnonorm, chanstoinclude, labels, xchan, labelclr );drawnow
        if 0
            figname = 'fig_nirshbo_pcmaps';  set( gcf, 'PaperPositionMode', 'auto');
            print(  '-dtiff', '-zbuffer', '-r600',  [figspath  figname  '.tiff'] );
        end
        plot6maps_horiz( 130, 'HbR', evecnr, lamnrnorm, chanstoinclude, labels, xchan, labelclr );drawnow
        if 0
            figname = 'fig_nirshbr_pcmaps';  set( gcf, 'PaperPositionMode', 'auto');
            print(  '-dtiff', '-zbuffer', '-r600',  [figspath  figname  '.tiff'] );
        end
    end
    %%
    if 0
        for ivecnum = 1:6,
            fh = figure(1234+ivecnum);clf;
            set(fh,'color','w');
            iverbose = 0; labelclr = 'k';
            rendermontage3( 1234, labels, xchan, iverbose, labelclr, chanstoinclude, evece(:,ivecnum), .1 );
            title(['Mode ' num2str(ivecnum) ': \lambda=' num2str(lamenorm(ivecnum)) ],'fontsize',14);
            for iii = 1:length(chanstoinclude), fprintf('[%s] %f\n',chanstoinclude{iii},evece(iii,ivecnum));end
        end
    end
    %% LAGGD CORR
    nlags = ceil( Fsnirs.*maxlagtime ); % max number of lags to compute (samples)
    ymin = -.3; ymax = .3; xmin = -maxlagtime; xmax = maxlagtime;
    fsz = 16;
    for jj = 2, % EEG
        for ii = 2, % NIRS
            xcorrelbar = [ ];
            xcorrelRbar = [ ];
            iend = 0;
            for kexp = 1:length(numsamples),
                ibeg = iend + 1;
                iend = iend + numsamples(kexp);
                picksamples = [ibeg:iend];
                aaae = aae(picksamples,jj);
                aaano = aano(picksamples,ii);
                aaanr = aanr(picksamples,ii);
                
                [ xcorrel, Lags ] = xcorr(  aano(:,ii), aae(:,jj),  nlags, 'coef' );
                [ xcorrelR, Lags ] = xcorr(   aanr(:,ii), aae(:,jj), nlags, 'coef' );
                xcorrelbar = [xcorrelbar xcorrel];
                xcorrelRbar = [xcorrelRbar xcorrelR];
            end
            fh = figure(4010);clf;hold on; set(fh,'color','w');
            plot(Lags./Fsnirs,xcorrelbar,'r','linewidth',3); plot(Lags./Fsnirs,xcorrelRbar,'b','linewidth',3);
            %title([ 'EEG-Hb Modes: ' num2str(jj) '-' num2str(ii) ', Minute ' num2str(kminutes) ]);
            titlestrn = ['EEG PC' num2str(jj) ', HbO PC' num2str(ii) ', HbR PC', num2str(ii)];
            title(titlestrn,'fontsize',fsz);
            axis([xmin xmax ymin ymax]); grid on;
            xlabel('Lag (s)','fontsize',fsz); ylabel('Correlation','fontsize',fsz); box on;
            set(gca,'fontsize',fsz);
            adfs=4;
            pause(2);
            
            if 0
                figname = ['fig_eeg_freq' num2str(ifreqband) '_eegpc' num2str(jj) '_hbpc' num2str(ii) '_'   ];
                set( gcf, 'PaperPositionMode', 'auto');
                print(  '-dtiff', '-zbuffer', '-r600',  [figspath  figname  '.tiff'] );
            end
        end
    end
    %% PC time course
    ai = aae(:,1);
    ttt = [1:length(aae(:,1))]./Fsnirs;
    fh = figure(5010);clf;hold on; set(fh,'color','w');
    plot(ttt,aae(:,1),ttt,aae(:,2)+1900);%,ttt,aae(:,3)+0,ttt,aae(:,4)+0); grid on;
    %axis([0 300 -1000 5000])
    
    for ii=1:16,
        for jj=ii:16,
            fprintf('%d %d %f\n',ii,jj,         sum( aae(:,ii).*aae(:,jj) ) );
        end
    end
    











end
