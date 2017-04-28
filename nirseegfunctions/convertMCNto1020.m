function chanlabels = convertMCNto1020(chanlabels)

for ii = 1:length(chanlabels),
    
    chanlabel = upper( char( chanlabels(ii) ) );
    
    if strcmp( chanlabel, 'T7' ),
        chanlabels{ii} = 'T3';
    end
    if strcmp( chanlabel, 'T8' ),
        chanlabels{ii} = 'T4';
    end
    if strcmp( chanlabel, 'P7' ),
        chanlabels{ii} = 'T5';
    end
    if strcmp( chanlabel, 'P8' ),
        chanlabels{ii} = 'T6';
    end
    
    
end



