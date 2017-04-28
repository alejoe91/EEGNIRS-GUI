function retval = getfieldval( filename, fieldname, delimiter ),

fid = fopen(filename);

tline = fgets(fid);
while ischar(tline)
    %disp(tline)
    tline = fgets(fid);
    charpos = strfind(tline, fieldname);
    %fprintf( '[%s] %d\n',tline,charpos);
    if ~isempty(charpos),
        delimpos = strfind(tline, delimiter);
        if ~isempty(delimpos),
            %fprintf('[%s]\n', tline(delimpos:end));
            retval = tline( delimpos+1:end );
            npos = strfind(retval, ';');
            retval(npos) = [];
            fclose(fid);
            return;
        end
    end
    
end

fclose(fid);






