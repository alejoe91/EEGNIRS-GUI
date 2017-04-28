function [eeg_chan_array, lab] = extract_eeg_array(configfile, eegchanlabels)

%this function reads the config file and returns the label, if used,
%corresponding to Source = isrc, Det = idet

% Read CONFIGURATION FILE

MAX_CHAN = 40;

fid = fopen(configfile);

tline = fgets(fid);
fprintf('EEG Configuration: %s',tline);

lab = cell(MAX_CHAN,1);

count = 1;

%Skip first line
    tline = fgets(fid);   
    
%load first line
    tline = fgets(fid);
    
while ischar(tline)         
    
    C = textscan(tline,'%s');
    lab(count) = C{1};
    count = count + 1;
    tline = fgets(fid);
    
end
   
lab = lab(1:count - 1);

%find the channels position in eegchanlabels and return the position array:

eeg_chan_array = zeros(length(lab),1);

found = [];

for i=1:length(lab)
    ind = find(strcmp(lab{i},eegchanlabels));
    if ~isempty(ind)
        eeg_chan_array(i) = ind;    
        found = [found , i];
    end
end

eeg_chan_array = eeg_chan_array(find(eeg_chan_array~=0));
lab = lab(found);

fclose(fid);

return