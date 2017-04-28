function eegpower = calc_eegpower_avgoverchans( chanstoinclude, flag_userelativeeegpower, ...
    numchaneeg,eegchanlabels,numfreqbands,xxeeg )

eegpower = cell( length(numfreqbands), 1 );

% freq bands: 1: 0-80,  2: 0-4,  3: 4-8,  4: 8-12,  5: 12-30,  6: 30-80,  7: 30-50,  8: 50-80

for ifreqband = 1:numfreqbands,
    eegdatabucket = [];
    for jj = 1:numchaneeg, % EEG
        chan2 = upper ( char( eegchanlabels(jj) ) );
        indxx = findstrincellarray( chanstoinclude, upper(chan2) );
        if indxx,
            eegdata = ( xxeeg(:,jj,ifreqband) );
            eegdatabucket = [eegdatabucket eegdata];
            fprintf('[%s]\n',chan2 );
        end
    end
    % avg over chans
    eegpower{ifreqband} = mean( eegdatabucket' );
end
% relative power
if flag_userelativeeegpower,
    eegtotalpower = eegpower{1};
    for ifreqband = 2:numfreqbands,
        eegdata = eegpower{ifreqband};
        eegdata = eegdata./eegtotalpower;
        eegpower{ifreqband} = eegdata;
    end
end