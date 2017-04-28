function chanlabel = get_chanlabel_from_srcdec( isrc, idet )

if isrc~=idet, error( 'UNEQUAL SRC AND DET NUMBER NOT SUPPORTED\n'); end

switch isrc,
    case 1, chanlabel = 'FP1';
    case 2, chanlabel = 'F3';
    case 3, chanlabel = 'F7';
    case 4, chanlabel = 'T3';
    case 5, chanlabel = 'C3';
    case 6, chanlabel = 'T5';
    case 7, chanlabel = 'P3';
    case 8, chanlabel = 'O1';
    case 9, chanlabel = 'FZ';
    case 10,chanlabel = 'CZ';
    case 11,chanlabel = 'PZ';
    case 12,chanlabel = 'FP2';
    case 13,chanlabel = 'F4';
    case 14,chanlabel = 'F8';
    case 15,chanlabel = 'C4';
    case 16,chanlabel = 'T4';
    case 17,chanlabel = 'P4';
    case 18, chanlabel = 'T6';
    case 19, chanlabel = 'O2';
    otherwise,      chanlabel = '***ERROR***';
end
return;
