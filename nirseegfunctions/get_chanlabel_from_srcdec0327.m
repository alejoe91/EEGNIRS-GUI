function chanlabel = get_chanlabel_from_srcdec( isrc, idet )

if isrc~=idet, error( 'UNEQUAL SRC AND DET NUMBER NOT SUPPORTED\n'); end

switch isrc,
    case 1, chanlabel = 'FP1';
    case 2, chanlabel = 'F3';
    case 3, chanlabel = 'F7';
    case 4, chanlabel = 'T7';
    case 5, chanlabel = 'P7';
    case 6, chanlabel = 'C3';
    case 7, chanlabel = 'P3';
    case 8, chanlabel = 'O1';
    case 9, chanlabel = 'FP2';
    case 10,chanlabel = 'FZ';
    case 11,chanlabel = 'CZ';
    case 12,chanlabel = 'PZ';
    case 13,chanlabel = 'O2';
    case 14,chanlabel = 'P4';
    case 15,chanlabel = 'C4';
    case 16,chanlabel = 'F4';
    case 17,chanlabel = 'F8';
    case 18, chanlabel = 'T8';
    case 19, chanlabel = 'P8';
    otherwise,      chanlabel = '***ERROR***';
end
return;
