function xxx = nirs_compressdata( mydat, chansrcdetindx, numsrc )
[nrows,ncols] = size( mydat );
numchan = length(chansrcdetindx);
xxx = zeros( nrows, numchan );

for ichan = 1:numchan,
    isrc = chansrcdetindx(ichan,1);
    idet = chansrcdetindx(ichan,2);
    icol = (isrc-1).*numsrc + idet;
    xxx(:,ichan) = mydat(:,icol);
end



