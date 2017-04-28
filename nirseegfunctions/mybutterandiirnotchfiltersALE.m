function xx = mybutterandiirnotchfiltersALE(xx, flags, order, Fs ,signal)
%
% xx contains a signal in each column
% flags = [ i j k l]
% i = 0: no highpass filtering, i = 0.23: highpass cutoff at 0.23 Hz
% j similarly for lowpass
% k similarly for notch, if k>0, then k is the frequency of notch
% l = 1: filtfilt, ie zero phase shift; 0: just filt
% order: butter
%

[nrow,ncol] = size(xx);
if nrow<ncol, fprintf('ERROR [%s]: each row of input data matrix must contain a signal time series\n',mfuncname); return; end

fhighpass = flags(1);
flowpass = flags(2);
fnotch = flags(3);

if strcmp(signal,'eeg')
    delf = .5; % notch half width;
else
    delf = .025; % notch half width;
end

if 0==fhighpass & 1000000<flowpass & 0==fnotch,
    return;
end

if length(flags)>3,
    if 0 < flags(4),
        zerophaseshift = 1;
    else,
        zerophaseshift = 0; 
    end
end

if 0 < fhighpass,
    [b,a] = butter( order, fhighpass ./(Fs./2), 'high' ); %    fvtool(b,a);
    %zerophaseshift+=1-->filtfilt-->filter with linear phase --> no distorsion
    %in harmonic content
    if 0==zerophaseshift,
        xx = filtfilt(b,a, xx );
    else
        xx = filter(b,a,xx);
    end
end

if 0 < flowpass,
    [b,a] = butter( order, flowpass ./(Fs./2), 'low' ); %    fvtool(b,a);
    if 0==zerophaseshift,
        xx = filtfilt(b,a, xx );
    else
        xx = filter(b,a,xx);
    end
end

if 0 < fnotch,
    %        wo = fnotch/(Fs/2);  bw = wo/35;        [b,a] = iirnotch(wo,bw); % fvtool(b,a);
    
    fnotch1 = (fnotch- delf)  ./(Fs./2); fnotch2 = (fnotch+ delf)  ./(Fs./2);
    [b,a] = butter( order, [fnotch1 fnotch2], 'stop' ); % fvtool(b,a);
    
    %{
[z p k] = butter(order, [fnotch1 fnotch2], 'stop'); % 10th order filter
[sos,g]=zp2sos(z,p,k); % Convert to 2nd order sections form
h=dfilt.df2sos(sos,g); % Create filter object
    %}
    
    if 0==zerophaseshift,
        xx = filtfilt(b,a, xx );
        %        xx = filter(h, xx );
    else
        xx = filter(b,a,xx);
    end
end
