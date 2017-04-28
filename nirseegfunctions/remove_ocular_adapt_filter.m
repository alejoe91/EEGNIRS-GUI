function [xc] = remove_ocular_adapt_filter( x, x1, x2, fs )

% Adaptive filetring

% "pathway" is the original folder, the data stored in, for example C:\Users\Desktop\Data\EEG
% "file" is the name of the file, for example Subject1.mat
% "x" is the original recorded signal containing 2 EOG channels at the 2 last columns,
% each channel will be placed in a column.

% "n" is total number of channels containing 2 EOG channels
% "fs" is sampling frequency
% "x1" is the signal from the first EOG channel
% "x2" is the signal from the second EOG channel
% "Clear" is the output without EOG channels
% "Clear2" is the output with EEG channels
%clear all
%close all
%path('pathway')
%load file
%fs= ;
[nrows,ncols] = size(x);
Clear=zeros( size(x) );

for i=1:ncols,
    d=x(:,i);            % Primary input
%    t = [0:length(d)-1]'/fs;
    % *** (USING TWO STAGES) *** %
%    N = length(d);
%    x1 = x(:,n-1);
%    x2 = x(:,n);
    % First filter
    mu = .0001;
    ha = adaptfilt.sd(32,mu);
    [y2,e2] = filter(ha,x2,d);
    % Second filter
    mu = .0001;
    ha = adaptfilt.sd(32,mu); % DSP toolbox
    [y3,e3] = filter(ha,x1,e2);
    
    Clear(:,i)=e3(:);
end
%Clear2=[Clear,x1,x2];
%xc=[Clear,x1,x2];
xc = Clear;
