function indx = findstrincellarray( cellarray, strn )

indx = 0;

N = length(cellarray);

for kk=1:N,
    
    sss = strtrim( char(cellarray{kk}) );   
        
    if strcmp( upper(sss), upper(strn) ),
        indx = kk;
        return;
    end
    

    
end


