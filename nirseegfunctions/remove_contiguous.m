function nevt = remove_contiguous( nevt0 )

indxremoveflag = [];

for ii = 2:length(nevt0),
    
    if nevt0(ii)==nevt0(ii-1)+1,
        indxremoveflag = [indxremoveflag; ii];
    end
    
end

nevt0(indxremoveflag) = [];

nevt = nevt0;


