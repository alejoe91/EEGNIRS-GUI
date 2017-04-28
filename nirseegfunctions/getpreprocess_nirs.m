function [nirschanlabel, xxnirsf ] = getpreprocess_nirs(filenirs,filterlofreq,filterhifreq, filternotch, Fsnirs)

% LOAD DATA
mydat = load( [ filenirs ]  );
% MAP NIRS SRC-DEC PAIR TO 10-20 CHAN LABEL
nirschanlabel = cell(19,1);
for ii = 1:19,
    nirschanlabel{ii} = get_chanlabel_from_srcdec( ii, ii );
    fprintf( '%d %d  %s\n', ii, ii, nirschanlabel{ii} );
end
% initialize the names and positions of predetermined full set of 128 channels
[ labels, xchan ] = inputchanlocations( 'positions_wholehead_flat.txt' );

% REDUCE DATA - OMIT UNUSED CHANS
%
% NOTES: 
% THERE IS A CORRECTION HERE FROM 20-20 TO 6-6
% ONLY 1-1, 2-2, 3-3, etc.
%
[nsamples, nchans0] = size( mydat );
xxnirs = zeros( nsamples, 19 );
for isrc = 1:20,
    for idet = 1:20,
        if isrc==idet, % we used 1-1, 2-2, etc.
            ichan = (isrc - 1).*20 + idet;
            if 20==isrc, % we used 20-20 for chan 6
                xxnirs( : , 6 ) = mydat( : , ichan );
            else, % we used 1-1 for chan1, 2-2 for chan2 etc.
                xxnirs( : , isrc ) = mydat( : , ichan );
            end
        end
    end
end
% FILTER
xxnirsf = mybutterandiirnotchfilters( xxnirs, [filterlofreq filterhifreq filternotch 1], 6, Fsnirs );
