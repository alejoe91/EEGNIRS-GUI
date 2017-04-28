function [genheadr1, sigheadr1, totnumrecs1, xx, Fs, numchan] = readEDFfile2ALE( file1, numrecordstoread, numchanmax );
% this is a wrapper around readEDFfile
% assumes all channels have same sample rate

[genheadr1, sigheadr1, totnumrecs1, data1] = readEDFfile( file1, numrecordstoread );

% sample freq
kk = 1;
nsamplesperrecord = str2num(sigheadr1(kk).Number_of_Samples);
durrecord = genheadr1.durationrec;
Fs = nsamplesperrecord./durrecord;

% num channels
numchan = genheadr1.nsig;

if numchanmax<numchan,
    numchan = numchanmax;
end

xx = zeros(length(data1{1}),numchan); 

% assign the data
for channum=1:numchan,
    xx(:,channum) = data1{channum}(1:length(data1{1}));
    %xx(:,channum) = data1{channum}(1:end);
    a = sigheadr1(channum);
    pmin = str2num(a.physical_minimum);
    pmax = str2num(a.physical_maximum);
    dmin =  str2num(a.digital_minimum);
    dmax =  str2num(a.digital_maximum);
    dp = pmax - pmin;
    dd = dmax - dmin;
    factr = dp./dd;
    pmid = .5.*(pmax+pmin);
    xx(:,channum) = pmid + factr .* xx(:,channum);
%   fprintf('\n%d  %f %f\n', channum, pmin, factr );    fprintf('%d %d %d %d\n', pmin, pmax, dmin, dmax);
end
