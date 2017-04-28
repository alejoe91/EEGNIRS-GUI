function         renderheaddots( chanlabels, labelfontsize, fontclr, xxx, yyy, zzz, titlestr, titlefontsz, clmap, dotsiz, flag_showaxis, valmax, valmin )

mapsz = length(clmap);
hold on;
offsetx = .6; offsety = .6;


axis( [ min(xxx) max(xxx) min(yyy) max(yyy)  ] );
for kk = 1:length(xxx),
    val = zzz(kk);
    
    if val<valmin,
        clrindx = 1;
    else,
        clrindx = 1 + floor(  (val - valmin)./(valmax-valmin) .* mapsz );
        clrindx = min( clrindx, mapsz );
    end
    mycolr = clmap( clrindx, : );
    
    plot( xxx(kk), yyy(kk), 'o','markerfacecolor',mycolr,'markeredgecolor',mycolr,'markersize', dotsiz );
    
end

if labelfontsize,
    for kk = 1:length(xxx),
        myx = xxx(kk) + offsetx;
        myy = yyy(kk) + offsety;
        text( myx, myy, chanlabels{kk}, 'fontsize', labelfontsize,'color', fontclr );
    end
end


if ~isempty(titlestr),
    title( titlestr,'fontsize',titlefontsz );
end

if 0==flag_showaxis,
    axis off;
end