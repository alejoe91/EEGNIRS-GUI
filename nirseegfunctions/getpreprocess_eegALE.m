function [chanlabels, chanlocs, xx ,xxfeeg, evteeg, Fs] = getpreprocess_eegALE( file1 ...
    , dtbeg ...
    , dtend ...
    , filterlofreq ...
    , filterhifreq ...
    , filternotch );

% READ DATA
numrecordstoread = []; % to read all records set numrecordstoread=[]
numchanmax = 26; %21;
[genheadr1, sigheadr1, totnumrecs1, xx, Fs, numchan] = readEDFfile2ALE( file1, numrecordstoread, numchanmax );
% CLIP BEGINNING AND END
nbeg = Fs.*dtbeg;
nend = Fs.*dtend;
xx = xx( nbeg:end-nend, : );

% FILTER
xxf = mybutterandiirnotchfiltersALE( xx, [filterlofreq filterhifreq filternotch 0], 6, Fs,'eeg' );
% initialize the names and positions of predetermined full set of 128 channels
[ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );
subsetindxs = [];
%chanlabels = cell(numchan,1);chanlocs = zeros(numchan,3);
%fprintf('Channels Used:\n');
count = 1;
for ii = 1:numchan,
    chan1 = sigheadr1(ii).label;
    [xxx1, iii1 , ifound1 ] = getchanlocation( chan1, labels, xchan );
    if 1 | ifound1,
        subsetindxs = [subsetindxs; ii];
        chanlabels{count} = chan1;
        chanlocs(count,:) = xxx1;
        count = count + 1;
        %      fprintf( '%d %s\n', ii, chan1);
    end
end

%chanlocs = [];
xxfeeg = xxf(:, subsetindxs);

% use nonfiltered for event channel
evteeg = xx(:,1);