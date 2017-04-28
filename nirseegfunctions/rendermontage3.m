function rendermontage3( fignum, labels, xchan, iverbose, labelclr, chanstoinclude, evec, cutoff, displaylabel ),

red = [1 0 0];
blu = [0 0 1];
blk = [0 0 0];
wht = [1 1 1];
map = colormap;
mapsz = length(map);

hold on; axis( [ -4 4 -4 5  ] );

offsetx = .2; offsety = .2;

for ii = 1:19,
    chan1 = get_chanlabel_from_srcdec( ii, ii );
    [ xxx, iii , ifound ] = getchanlocation( chan1, labels, xchan );
    if iverbose, fprintf( '%d %s %d \n', ii, chan1, ifound ); end
    if ifound,
        
        plot( xxx(1), xxx(2), 'ro' );
        
        % check if this chan makes the cutoff
        ind = findstrincellarray( chanstoinclude, chan1 );
        if ~isempty(ind) & 0<ind,
            val = evec(ind);
            if cutoff<abs(val),
                fi = (val-(-1))./2;
                %                mycolr = fi.*(red)./1 + (1-fi).*(blu)./1;
                mycolr = map( 1+floor( fi.*mapsz ), : );
                plot( xxx(1), xxx(2), 'o','markerfacecolor',mycolr,'markeredgecolor',mycolr,'markersize',14 );
            end
        end
        if displaylabel,
            myx = xxx(1) + offsetx; myy = xxx(2) + offsety;
            text( myx, myy,  chan1, 'color', labelclr );
        end
    end
end
axis off