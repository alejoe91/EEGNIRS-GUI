function rendermontage( fignum, sigheadr1, labels, xchan, iverbose, labelclr ),

hold on; axis( [ -4 4 -4 5  ] );

offsetx = .2; offsety = .2;

for ii = 1:length( sigheadr1 ),
    chan1 = sigheadr1(ii).label;
    [ xxx, iii , ifound ] = getchanlocation( chan1, labels, xchan );
    if iverbose, fprintf( '%d %s %d \n', ii, chan1, ifound ); end
    if ifound,
        plot( xxx(1), xxx(2), 'ro' );
        myx = xxx(1) + offsetx; myy = xxx(2) + offsety;
        text( myx, myy,  chan1, 'color', labelclr );
    end
end
axis off