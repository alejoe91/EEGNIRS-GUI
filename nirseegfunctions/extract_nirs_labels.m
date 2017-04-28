function [src_det_pair, lab] = extract_nirs_labels(configfile)

%this function reads the config file and returns the label, if used,
%corresponding to Source = isrc, Det = idet

% Read CONFIGURATION FILE

MAX_CHAN = 40;

fid = fopen(configfile);

tline = fgets(fid);
fprintf('Configuration: %s',tline);

src_det_pair = zeros(MAX_CHAN,2);

lab = cell(MAX_CHAN,1);

count = 1;

%Skip first line
    tline = fgets(fid);   
    
%load first line
    tline = fgets(fid);
    
while ischar(tline)         
    
    C = textscan(tline,'%d\t%d\t%s');
    src_det_pair(count,1) = C{1} ;
    src_det_pair(count,2) = C{2};
    lab(count) = C{3};
    count = count + 1;
    tline = fgets(fid);
    
end
   
src_det_pair = src_det_pair(1:count-1,:);
lab = lab(1:count - 1);

fclose(fid);

return