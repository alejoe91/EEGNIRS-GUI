function xxnirs = fnnirs_reduce_data_2( numchan, mydat )
% picks out only those chans where src# == det#

[nsamples, dumnum] = size( mydat );

xxnirs = zeros( nsamples, numchan);

for isrc = 1:numchan,
    for idet = 1:numchan,
        if isrc==idet, % we used 1-1, 2-2, etc.
            ichan = (isrc - 1).*numchan + idet;
            % we used 1-1 for chan1, 2-2 for chan2 etc.
            xxnirs( : , isrc ) = mydat( : , ichan );
        end
    end
end