function rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr ),

%rendermontage( fignum_meanadjmat, sigheadr1, labels, xchan, iverbose, labelclr );

figure(fignum);clf; hold on; axis( [ -4 4 -4 4  ] );
gry = [.6 .6 .6];
offsetx = .2; offsety = .2;

for ii = 1:numchan,
    chan1 = sigheadr1(ii).label;
    [ xxx, iii , ifound ] = getchanlocation( chan1, labels, xchan );
    fprintf( '%d %s %d \n', ii, chan1, ifound );
    if ifound,
        plot( xxx(1), xxx(2), 'ro' );
        myx = xxx(1) + offsetx; myy = xxx(2) + offsety;
        text( myx, myy,  chan1, 'color', gry );
    end
end
