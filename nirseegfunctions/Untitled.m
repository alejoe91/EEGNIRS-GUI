clear;
%% input labels and their positions 

fid = fopen('positions.txt');
tline = fgetl(fid);
count = 1;
while ischar(tline) & ~isempty(tline),
    C = textscan(tline,'%s %f %f %f');    
    label{count} = char(C{1}); 
    xx(count,1) = C{2};
    xx(count,2) = C{3};
    xx(count,3) = C{4};
    count = count + 1;
    tline = fgetl(fid)
end
fclose(fid)

%%

