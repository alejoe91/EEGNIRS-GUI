
function xxnirs = fnnirs_reduce_date( mydat )

[nsamples, nchans0] = size( mydat );

numchans = sqrt(nchans0);

xxnirs = zeros( nsamples, numchans );

for isrc = 1:numchans,
    for idet = 1:numchans,
        if isrc==idet, % we used 1-1, 2-2, etc.
            ichan = (isrc - 1).*numchans + idet;
            % we used 1-1 for chan1, 2-2 for chan2 etc.
            xxnirs( : , isrc ) = mydat( : , ichan );
        end
    end
end