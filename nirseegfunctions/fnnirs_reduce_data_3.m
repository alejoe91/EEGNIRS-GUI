function     xx = fnnirs_reduce_data_3( icycle, numsrc, numdet, mydat )

if 1122==icycle | 1123==icycle | 1132==icycle | 1133==icycle | 1134==icycle| 1135==icycle...
        | 1142==icycle| 1143==icycle| 1144==icycle| 1145==icycle...
         | 1152==icycle| 1153==icycle| 1154==icycle| 1155==icycle,
    
    countchan = 1;
    for isrc = 1:numsrc,
        for idet = 1:numdet,            
            ichan = (isrc - 1).*numdet + idet;            
            tmpstrng = get_chanlabel_from_srcdetpair_multi( icycle, isrc, idet );            
            if isempty(findstr(tmpstrng, 'UNUSED')),
                fprintf('isrc,idet %d,%d   ichan %d   countchan %d   chanlabel [%s]\n',isrc,idet,ichan,countchan, tmpstrng);  
%mydat contains all the possible numsrc x numdet combination
                xx( : , countchan ) = mydat( : , ichan );
                countchan = countchan + 1;
            end            
        end
    end
    return;
end
        
    
if 1102==icycle | 1103==icycle...
        | 1112==icycle | 1113==icycle
    
    % experiments where src index equals det index
    xx = fnnirs_reduce_data_2( numchan, mydat );
    return;
end







