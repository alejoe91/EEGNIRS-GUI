function [labels, xx] = inputchanlocations( fname );

fid = fopen(  fname  );
tline = fgetl(fid);
count = 1;
while ischar(tline) & ~isempty(tline),
    C = textscan(tline,'%s %f %f %f');
    labels{count} = char(C{1});
    xx(count,1) = C{2};
    xx(count,2) = C{3};
    xx(count,3) = C{4};
    count = count + 1;
    tline = fgetl(fid);
end
fclose(fid);


% insert additional positions on right side
[xx, labels] = insertadditionalpositions( xx, labels );

% insert left side by taking the negative of x axis
[xx, labels] = insertleftpositions( xx, labels );
