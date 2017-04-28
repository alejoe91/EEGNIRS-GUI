function [ xxx, iii, ifound ] = getchanlocation( chanstrinput, labels, xchan )

chanstrinput = upper( chanstrinput );

for kk = 1:length(labels),

    ccc = upper( char( labels(kk) ) );

    if strcmp( ccc, chanstrinput ),
        iii = kk;
        xxx = xchan( kk, : );
        ifound = 1;
        return
    end

end

xxx = 999.*ones(1,3);
iii = 999;
ifound = 0;


