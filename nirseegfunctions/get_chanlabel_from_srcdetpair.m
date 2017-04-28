function chanlabel = get_chanlabel_from_srcdetpair( isrc, idet )

if isrc~=idet, error( 'UNEQUAL SRC AND DET NUMBER NOT SUPPORTED\n'); end

switch isrc,
    case 1, chanlabel = 'FP1';
    case 2, chanlabel = 'F3';
    case 3, chanlabel = 'F7';
    case 4, chanlabel = 'T7';
    case 5, chanlabel = 'T5';
    case 6, chanlabel = 'C3';
    case 7, chanlabel = 'FZ';
    case 8, chanlabel = 'O1';
    case 9, chanlabel = 'FP2';
    case 10,chanlabel = 'CZ';
    case 11,chanlabel = 'F4';
    case 12,chanlabel = 'C4';
    case 13,chanlabel = 'O2';
    case 14,chanlabel = 'F8';
    case 15,chanlabel = 'T8';
    case 16,chanlabel = 'T6';
    otherwise, chanlabel = '***ERROR***';
end
return;

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
    otherwise, chanlabel = '***ERROR***';
end
return;

%{
16 channels list - 
Fp1- Ch1 ( S1-D1
F3 – Ch2 
F7 – Ch3
T7 – Ch4
T5- Ch5
C3 - Ch6
Fz- Ch7
O1- Ch8
Fp2 – Ch9
Cz – Ch10
F4- Ch11
C4 –Ch12
O2 –Ch`13
F8-Ch14
T8 – Ch15
T6 – Ch 16
%}

%{
ordered as they appear in the corr matrix
FP1 FP2 F7 F8 F3 F4 FZ T7 T8 C3 C4 CZ P7 P8 P3 P4 PZ O1 O2
1   9   3  17 2  16 10 4  18 6  15 11 5  19 7  14 12 8  13
FP1 
FP2
F7
F8
F3
F4
FZ
T7
T8
C3
C4
CZ
P7
P8
P3
P4
PZ
O1
O2
%}