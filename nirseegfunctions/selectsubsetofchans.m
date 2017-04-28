function  xxret = selectsubsetofchans( xxe, chanlabels, chanstoinclude );
ind = [];
for iii = 1:length(chanstoinclude),
    mychan = upper( chanstoinclude{iii} );
    ind = [ind; findstrincellarray( chanlabels, mychan ) ];
end
xxret = xxe(:,ind);