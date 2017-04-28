function chanlabel = get_chanlabel_from_srcdetpair_multi( iflag, isrc, idet )

if 1102==iflag | 1103==iflag...
        | 1112==iflag | 1113==iflag
    
    if isrc~=idet, error( 'FOR THIS EXPERIMENT SRC-DET index NUMBERS MUST BE EQUAL\n'); end
    
    switch isrc,
        case 1, chanlabel = 'FP1';
        case 2, chanlabel = 'F3';
        case 3, chanlabel = 'F7';
        case 4, chanlabel = 'T7';
        case 5, chanlabel = 'C3';
        case 6, chanlabel = 'P7';
        case 7, chanlabel = 'O1';
        case 8, chanlabel = 'O2';
        case 9, chanlabel = 'P3';
        case 10,chanlabel = 'PZ';
        case 11,chanlabel = 'P4';
        case 12,chanlabel = 'SS'; %% superficial src-det pair
        case 17,chanlabel = 'FP2';
        case 18,chanlabel = 'FZ';
        case 19,chanlabel = 'CZ';
        case 20,chanlabel = 'F4';
        case 21,chanlabel = 'F8';
        case 22,chanlabel = 'T8';
        case 23,chanlabel = 'C4';
        case 24,chanlabel = 'T6';
        otherwise, chanlabel = '***ERROR***';
    end
    
elseif 1122==iflag | 1123==iflag | 1132==iflag | 1133==iflag | 1134==iflag...
        | 1135==iflag | 1142==iflag | 1143==iflag | 1144==iflag| 1145==iflag...
        | 1152==iflag | 1153==iflag | 1154==iflag| 1155==iflag, 
    
    % Unequal src-det index numbers are superficial chans
    
    if 1==isrc & 1==idet, chanlabel = 'FP1';
    elseif 2==isrc & 2==idet, chanlabel = 'F3';
    elseif 3==isrc & 3==idet, chanlabel = 'F7';
    elseif 4==isrc & 4==idet, chanlabel = 'T3';
    elseif 5==isrc & 5==idet, chanlabel = 'C3';
    elseif 6==isrc & 6==idet, chanlabel = 'T5';
    elseif 7==isrc & 7==idet, chanlabel = 'P3';
    elseif 8==isrc & 8==idet, chanlabel = 'O1';
    elseif 9==isrc & 9==idet, chanlabel = 'FP2';
    elseif 10==isrc & 10==idet,chanlabel = 'F4';
    elseif 11==isrc & 11==idet,chanlabel = 'F8';
    elseif 12==isrc & 12==idet,chanlabel = 'T4';
    elseif 13==isrc & 13==idet,chanlabel = 'C4';
    elseif 14==isrc & 14==idet,chanlabel = 'T6';
    elseif 15==isrc & 15==idet,chanlabel = 'P4';
    elseif 16==isrc & 16==idet,chanlabel = 'O2';
    elseif 17==isrc & 17==idet,chanlabel = 'FZ';
    elseif 18==isrc & 18==idet,chanlabel = 'CZ';
    elseif 19==isrc & 19==idet,chanlabel = 'PZ';
    elseif 20==isrc & 20==idet,chanlabel = 'SSUPER';
    elseif 3==isrc & 23==idet, chanlabel = 'F7SUPER';
    elseif 4==isrc & 24==idet, chanlabel = 'T3SUPER';
    elseif 11==isrc & 21==idet, chanlabel = 'F8SUPER';
    elseif 12==isrc & 22==idet, chanlabel = 'T4SUPER';
    else, chanlabel = 'UNUSED';
    end
    
end

return;
