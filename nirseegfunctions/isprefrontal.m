function  retval = isprefrontal( chanlabel )

chanlabel = upper(chanlabel);

switch chanlabel,
    case 'FP1', retval = true;
    case 'F7', retval = true;
    case 'F3', retval = true;
    case 'FZ', retval = true;
    case 'FP2', retval = true;
    case 'F4', retval = true;
    case 'F8', retval = true;

    otherwise,  retval = false;
end