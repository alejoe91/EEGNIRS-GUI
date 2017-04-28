clear;
%% set folder for data
parentfolder = 'F:\LabData\LabData\';
subfolder = 'AhmetOmurtag0801\EMG\';
fnam = 'testdata.mat';
%% choose experiment
experimentnumber = 2;
%{
1: Resting
2: VFT pronouncing
3: VFT writing
4: Artifacts: jaw clench, rest, smile, rest, brow knit, rest; repeat
%}
%% load data
loadstring = ['load '  parentfolder subfolder fnam ];
fprintf('[%s]\n', loadstring);
eval( loadstring );
%% assign data
if 1==experimentnumber,
    myexperim = test01;
    experimentstring = ' Resting ';
elseif 2==experimentnumber,
    myexperim = test02;
    experimentstring = ' VFT pronouncing ';
elseif 3==experimentnumber,
    myexperim = test03;
    experimentstring = ' VFT writing ';
elseif 4==experimentnumber,
    myexperim = test04;
    experimentstring = ' Artifacts ';
else,
    error('Choose 1-4');
end
%% sample rate
Fs = myexperim.fs;
%% pick up first event
eventsig=myexperim.data{137};
thresh = -5;
indx = remove_contiguous(  find( thresh>eventsig )  );
%figure(5454);clf; plot(eventsig); hold on; plot(indx, eventsig(indx),'ro');
indxfirstevent = indx(1);
%% loop through chans
for ii = 1:4,
    %{
    Channel 1 – Temporalis 1
    Channel 2 – Temporalis 2
    Channel 3 - Masseter
    Channel 4 – Mentalis
    %}
    if 1==ii, strng = ' Temporalis 1 ';
    elseif 2==ii, strng = ' Temporalis 2 ';
    elseif 3==ii, strng = ' Masseter ';
    elseif 4==ii, strng = ' Mentalis ';
    else, strng = 'ERROR'; end
    
    aa = myexperim.data{ii};
    % clip at first event
    aa = aa(indxfirstevent:end);
    
    mytime = [1:length(aa)]./Fs;
    
    %% plot raw data
    figure(5456);clf;
    xmin = mytime(1); xmax = mytime(end);
    % xmin = 0; xmax = 870;
    ymin = min(aa); ymax = max(aa);
    plot(mytime,mysmooth(aa,200));
    axis([xmin xmax ymin ymax]);
    title(['EMG ' experimentstring strng]);
    
    %% power spectrogram
    window = ceil(Fs.*.8); % window size num samples
    noverlap = ceil(window./2);
    F = [0:1:800];
    fprintf('Calculating spectrogram...\n');
    [S,F,T,P] = spectrogram(aa, window, noverlap, F, Fs);
    % S: Fourier coeffs
    % F: freqs
    % T: times
    % P: PSD   
    %% gamma
    indx = find( 30<F & F<80 );
    Ptmp = P(indx,:);
    pp = sum(Ptmp,1);
    figure(7474); clf;
    plot( T, pp );
    title(['EMG gamma power (30-80Hz) ' experimentstring strng]);    
    drawnow;
    %% plot full spectrogram
     figure(6464);clf;
    surf(T,F,10.*log10(P),'edgecolor','none');
    view(0,90);
    xlabel('Time (s)'); ylabel('Hz');
    title(['EMG spectrogram ' experimentstring strng]);
    drawnow;
    %%
    retval = input('Press <Enter> to continue...');
        
end


