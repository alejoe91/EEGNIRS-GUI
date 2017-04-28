function [xx, label] = insertadditionalpositions( xx, label );

count = length(label) + 1;
label{count} = '74';xx(count,:) = ( xx(findstrincellarray( label, 'AFz' ),:) + xx(findstrincellarray( label, 'Fp2' ),:) ) ./ 2;count = count + 1;
label{count} = '77';xx(count,:) = ( xx(findstrincellarray( label, 'Fz' ),:) + xx(findstrincellarray( label, 'AF4' ),:) ) ./ 2;count = count + 1;
label{count} = '78';xx(count,:) = ( xx(findstrincellarray( label, 'F4' ),:) + xx(findstrincellarray( label, 'AF8' ),:) ) ./ 2;count = count + 1;
label{count} = '84';xx(count,:) = ( xx(findstrincellarray( label, 'FCz' ),:) + xx(findstrincellarray( label, 'F2' ),:) ) ./ 2;count = count + 1;
label{count} = '85';xx(count,:) = ( xx(findstrincellarray( label, 'FC2' ),:) + xx(findstrincellarray( label, 'F4' ),:) ) ./ 2;count = count + 1;
label{count} = '86';xx(count,:) = ( xx(findstrincellarray( label, 'FC4' ),:) + xx(findstrincellarray( label, 'F6' ),:) ) ./ 2;count = count + 1;
label{count} = '87';xx(count,:) = ( xx(findstrincellarray( label, 'FC6' ),:) + xx(findstrincellarray( label, 'F8' ),:) ) ./ 2;count = count + 1;
label{count} = '88';xx(count,:) = ( xx(findstrincellarray( label, 'FT8' ),:) + xx(findstrincellarray( label, 'F10' ),:) ) ./ 2;count = count + 1;
label{count} = '94';xx(count,:) = ( xx(findstrincellarray( label, 'Cz' ),:) + xx(findstrincellarray( label, 'FC2' ),:) ) ./ 2;count = count + 1;
label{count} = '95';xx(count,:) = ( xx(findstrincellarray( label, 'C2' ),:) + xx(findstrincellarray( label, 'FC4' ),:) ) ./ 2;count = count + 1;
label{count} = '96';xx(count,:) = ( xx(findstrincellarray( label, 'C4' ),:) + xx(findstrincellarray( label, 'FC6' ),:) ) ./ 2;count = count + 1;
label{count} = '97';xx(count,:) = ( xx(findstrincellarray( label, 'C6' ),:) + xx(findstrincellarray( label, 'FT8' ),:) ) ./ 2;count = count + 1;
label{count} = '98';xx(count,:) = ( xx(findstrincellarray( label, 'T8' ),:) + xx(findstrincellarray( label, 'FT10' ),:) ) ./ 2;count = count + 1;

label{count} = '103';xx(count,:) = ( xx(findstrincellarray( label, 'CPz' ),:) + xx(findstrincellarray( label, 'C2' ),:) ) ./ 2;count = count + 1;
label{count} = '104';xx(count,:) = ( xx(findstrincellarray( label, 'CP2' ),:) + xx(findstrincellarray( label, 'C4' ),:) ) ./ 2;count = count + 1;
label{count} = '105';xx(count,:) = ( xx(findstrincellarray( label, 'CP4' ),:) + xx(findstrincellarray( label, 'C6' ),:) ) ./ 2;count = count + 1;
label{count} = '106';xx(count,:) = ( xx(findstrincellarray( label, 'CP6' ),:) + xx(findstrincellarray( label, 'T8' ),:) ) ./ 2;count = count + 1;
label{count} = '112';xx(count,:) = ( xx(findstrincellarray( label, 'Pz' ),:) + xx(findstrincellarray( label, 'CP2' ),:) ) ./ 2;count = count + 1;
label{count} = '113';xx(count,:) = ( xx(findstrincellarray( label, 'P2' ),:) + xx(findstrincellarray( label, 'CP4' ),:) ) ./ 2;count = count + 1;
label{count} = '114';xx(count,:) = ( xx(findstrincellarray( label, 'P4' ),:) + xx(findstrincellarray( label, 'CP6' ),:) ) ./ 2;count = count + 1;
label{count} = '115';xx(count,:) = ( xx(findstrincellarray( label, 'P6' ),:) + xx(findstrincellarray( label, 'TP8' ),:) ) ./ 2;count = count + 1;
label{count} = '116';xx(count,:) = ( xx(findstrincellarray( label, 'P8' ),:) + xx(findstrincellarray( label, 'TP10' ),:) ) ./ 2;count = count + 1;
label{count} = '120';xx(count,:) = ( xx(findstrincellarray( label, 'POz' ),:) + xx(findstrincellarray( label, 'P2' ),:) ) ./ 2;count = count + 1;
label{count} = '121';xx(count,:) = ( xx(findstrincellarray( label, 'P4' ),:) + xx(findstrincellarray( label, 'PO8' ),:) ) ./ 2;count = count + 1;
label{count} = '125';xx(count,:) = ( xx(findstrincellarray( label, 'POz' ),:) + xx(findstrincellarray( label, 'O2' ),:) ) ./ 2;count = count + 1;
label{count} = '128';xx(count,:) = ( xx(findstrincellarray( label, 'Oz' ),:) + xx(findstrincellarray( label, 'O10' ),:) ) ./ 2;count = count + 1;
label{count} = '126';xx(count,:) = ( xx(findstrincellarray( label, 'O2' ),:) + xx(findstrincellarray( label, 'PO1' ),:) ) ./ 2;count = count + 1;
label{count} = '122';xx(count,:) = ( xx(findstrincellarray( label, 'PO8' ),:) + xx(findstrincellarray( label, 'P10' ),:) ) ./ 2;count = count + 1;

