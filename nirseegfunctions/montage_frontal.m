clear;
%% params
maxsepar = 2; % cutoff src-det distance for nirs channel
Ntrial = 1e2;
Maxtrial = 1e10;
labeloffset = .15;
%% max numbers of sensors
nsmax = 24; % source 1
ndmax = 24; % detector 2
nemax = 20; % eeg 3
%%
gry = [.6 .6 .6];
eeggry = [.4 .0 .4];
nirschangry = [.6 .6 .6];
srcolor = 'y';
detcolor = 'g';
eegcolor = eeggry;
%% input labels and their positions
fid = fopen('positions_prefrontal_right_74_flat.txt');
tline = fgetl(fid);
count = 1;
while ischar(tline) & ~isempty(tline),
    C = textscan(tline,'%s %f %f %f');
    label{count} = char(C{1});
    xx(count,1) = C{2};
    xx(count,2) = C{3};
    xx(count,3) = C{4};
    count = count + 1;
    tline = fgetl(fid);
end
fclose(fid);

%% additional positions
count = length(label) + 1;

label{count} = '74';
xx(count,:) = ( xx(1,:) + xx(17,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '77';
xx(count,:) = ( xx(9,:) + xx(15,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '78';
xx(count,:) = ( xx(3,:) + xx(7,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '84';
xx(count,:) = ( xx(12,:) + xx(22,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '85';
xx(count,:) = ( xx(18,:) + xx(3,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '86';
xx(count,:) = ( xx(19,:) + xx(2,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '87';
xx(count,:) = ( xx(20,:) + xx(6,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '88';
xx(count,:) = ( xx(21,:) + xx(11,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '94';
xx(count,:) = ( xx(10,:) + xx(18,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '95';
xx(count,:) = ( xx(14,:) + xx(19,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '96';
xx(count,:) = ( xx(4,:) + xx(20,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '97';
xx(count,:) = ( xx(5,:) + xx(21,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

label{count} = '98';
xx(count,:) = ( xx(13,:) + xx(8,:) ) ./ 2;
%plot( xx(count, 1), xx(count, 2), 'bo' );   text( xx(count, 1)+labeloffset, xx(count, 2), label( count ), 'color', gry );
count = count + 1;

%% 2d plot
if 1,
    figure(1010);clf; hold on;
    labeloffset = .2;
    for kk = 1:length(label),
        plot( xx(kk, 1), xx(kk, 2), 'ro' );
        text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', gry );
    end
    axis([-1 6 -1 6])
end
%% 3d plot
if 0,
    figure(1020);clf; hold on;
    for kk = 1:length(label),
        plot3( xx(kk, 1), xx(kk, 2), xx(kk, 3), 'ro' );
        text( xx(kk, 1)+labeloffset, xx(kk, 2), xx(kk, 3), label( kk ), 'color', gry );
    end
    axis([-1 14 -1 14 -1 14])
end

%% generate random montage
numpos = length(label);
itrial = 1;
tottrial = 1;
mysepareeg = 10000;
myseparnirschan = 10000;
myseparnirssrc = 10000;
myseparnirsdet = 10000;
mysepar = 10000;
countgood = 1;

while itrial <= Ntrial & tottrial <= Maxtrial,
    rrr = ceil(3.*rand( numpos,1));
    tottrial = tottrial + 1;

    % check if sensor numbers are under max
    indx = find( 1==rrr );
    if isempty(indx), continue; end;
    ns = length( indx );

    indx = find( 2==rrr );
    if isempty(indx), continue; end;
    nd = length( indx );

    indx = find( 3==rrr );
    if isempty(indx), continue; end;
    ne = length( indx );

    % nboundary = 9; ninterior = 26; % for the prefrontal region of choice
    factor = 1+26./35;
    %fprintf( 'ns, nd, ne %d %d %d\n',ceil(ns.*factor) ,ceil(nd.*factor),ceil( ne.*factor));

    if ceil(ns.*factor)>nsmax | ceil(nd.*factor)>ndmax | ceil(ne.*factor)>nemax, continue; end

    %% nirs channels
    if itrial>1,
        clear cc ccsrc ccdet
    end
    countchan = 1;
    for kk = 1: numpos,
        if 1==rrr(kk), % source
            for nn = 1 : numpos,
                if 2==rrr(nn), % detector
                    separ = sqrt(  ( xx(kk,1)-xx(nn,1)  ).^2 + ( xx(kk,2)-xx(nn,2)  ).^2  );
                    if maxsepar >= separ, % this is a channel
                        cc(countchan,1) = ( xx(kk,1)+xx(nn,1) )./2;
                        cc(countchan,2) = ( xx(kk,2)+xx(nn,2) )./2;
                        ccsrc(countchan,1) = xx(kk,1);
                        ccsrc(countchan,2) = xx(kk,2);
                        ccdet(countchan,1) = xx(nn,1);
                        ccdet(countchan,2) = xx(nn,2);
                        countchan = countchan + 1;
                    end
                end
            end
        end
    end
    numnirschan = countchan -1;
    %% plot
    if 0
        fh = figure(110); clf; set(fh, 'color', 'white');  hold on;
        %% plot nirs chans
        for kk = 1: numnirschan,
            pnirschan = plot( cc(kk, 1), cc(kk, 2), 'markersize',10,'marker','s','markerfacecolor', [nirschangry],'markeredgecolor',[nirschangry]);
            plot( [ccsrc(kk,1) ccdet(kk,1)],  [ccsrc(kk,2) ccdet(kk,2)] ,'color',  [nirschangry] );
        end
        %% plot src, det, eegs
        for kk = 1: numpos,
            if 1==rrr(kk), colr = srcolor;
            elseif 2==rrr(kk), colr = detcolor;
            elseif 3==rrr(kk), colr = eegcolor;
            else, error('should not be here'); end;
            if 1==rrr(kk),
                psrc = plot( xx(kk, 1), xx(kk, 2), 'markersize',6,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 2==rrr(kk),
                pdet = plot( xx(kk, 1), xx(kk, 2), 'markersize',6,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 3==rrr(kk), % eeg
                peeg = plot( xx(kk, 1), xx(kk, 2), 'markersize',10,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
            end
            text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', gry, 'fontsize',8 );
        end
        axis([-.2 6 -.2 5])
        box off
        axis off
        legend( [pnirschan psrc pdet peeg], 'NIRS Channel', 'Src', 'Det','EEG' );
    end
    %% extract eeg sensors
    indx = find( 3==rrr );
    xxeeg = xx( indx, : );
    labeleeg = label( indx );

    %% determine distance to nearest eeg sensor
    for kk = 1:numpos,
        distmineeg = 10000;
        for nn = 1:ne,
            dist = sqrt( (xx(kk,1)-xxeeg(nn,1)).^2 + (xx(kk,2)-xxeeg(nn,2) ).^2 );
            if dist < distmineeg,
                distmineeg = dist;
            end
        end
        disteeg(kk) = distmineeg;
    end
    % largest distance (from any position) to the nearest eeg sensor
    distmaxeeg = max( disteeg );
    %% determine distance to nearest nirs chan
    for kk = 1:numpos,
        distmin = 10000;
        for nn = 1:length( cc ),
            dist = sqrt( (xx(kk,1)-cc(nn,1)).^2 + (xx(kk,2)-cc(nn,2) ).^2 );
            if dist < distmin,
                distmin = dist;
            end
        end
        distnirs(kk) = distmin;
    end
    % largest distance (from any position) to the nearest nirs chan
    distmaxnirschan = max( distnirs );
    %% determine distance to nearest nirs source
    for kk = 1:numpos,
        distmin = 10000;
        for nn = 1:numpos,
            if 1==rrr(nn) & kk~=nn,
                dist = sqrt( (xx(kk,1)-xx(nn,1)).^2 + (xx(kk,2)-xx(nn,2) ).^2 );
                if dist < distmin,
                    distmin = dist;
                end
            end
        end
        distnirs(kk) = distmin;
    end
    % largest distance (from any position) to the nearest nirs chan
    distmaxnirssrc = max( distnirs );
    %% determine distance to nearest nirs det
    for kk = 1:numpos,
        distmin = 10000;
        for nn = 1:numpos,
            if 2==rrr(nn) & kk~=nn,
                dist = sqrt( (xx(kk,1)-xx(nn,1)).^2 + (xx(kk,2)-xx(nn,2) ).^2 );
                if dist < distmin,
                    distmin = dist;
                end
            end
        end
        distnirs(kk) = distmin;
    end
    % largest distance (from any position) to the nearest nirs chan
    distmaxnirsdet = max( distnirs );
    %%

    % crietrion
    %    ddd = sqrt(distmaxeeg.^2 + distmaxnirschan.^2);
    %    ddd = sqrt(distmaxeeg.^2 + distmaxnirssrc.^2 + distmaxnirsdet.^2);
    myvec = [distmaxeeg distmaxnirschan distmaxnirssrc distmaxnirsdet];
    ddd = sqrt( mean(myvec.^2) );
    % update minima
    %    if distmaxeeg<mysepareeg & distmaxnirs<myseparnirs,
    if ddd<mysepar,

        mysepar = ddd;
        mysepareeg = distmaxeeg;
        myseparnirschan = distmaxnirschan; myseparnirssrc = distmaxnirssrc; myseparnirsdet = distmaxnirsdet;
        rrrlist{countgood} = rrr;
        cclist{countgood} = cc;
        ccsrclist{countgood} = ccsrc;
        ccdetlist{countgood} = ccdet;
        myseparlist(countgood) = mysepar;
        mysepareeglist(countgood) = mysepareeg;
        myseparnirschanlist(countgood) = myseparnirschan;
        myseparnirssrclist(countgood) = myseparnirssrc;
        myseparnirsdetlist(countgood) = myseparnirsdet;
        countgood = countgood + 1;

        %% good plot
        fh = figure(220); clf; set(fh, 'color', 'white');  hold on;
        %% plot nirs chans
        for kk = 1: length(cc),
            pnirschan = plot( cc(kk, 1), cc(kk, 2), 'markersize',10,'marker','s','markerfacecolor', [nirschangry],'markeredgecolor',[nirschangry]);
            plot( [ccsrc(kk,1) ccdet(kk,1)],  [ccsrc(kk,2) ccdet(kk,2)] ,'color',  [nirschangry] );
        end
        %% plot eeg
        for kk = 1: numpos,
            if 1==rrr(kk), colr = srcolor;
            elseif 2==rrr(kk), colr = detcolor;
            elseif 3==rrr(kk), colr = eegcolor;
            else, error('should not be here'); end;
            if 1==rrr(kk),
                psrc = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 2==rrr(kk),
                pdet = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 3==rrr(kk), % eeg
                peeg = plot( xx(kk, 1), xx(kk, 2), 'markersize',8,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
            end
            text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', gry, 'fontsize',8 );
        end
        axis([-.2 6 -.2 5])
        box off
        axis off
        legend( [pnirschan psrc pdet peeg], 'NIRS Channel', 'Src', 'Det','EEG', 'fontsize',9 );
        legend boxoff;
        titlestr = [ ...
            num2str( itrial  ) ' of ' ...
            num2str( Ntrial ) ':  \Delta=' ...            
            num2str( mysepar ) '  \Delta_{E}=' ...
            num2str( mysepareeg  ) '   \Delta_{NC}=' ...
            num2str( myseparnirschan )  '   \Delta_{NS}=' ...
            num2str( myseparnirssrc )  '   \Delta_{ND}=' ...
            num2str( myseparnirsdet ) ];
        title( titlestr, 'fontsize',9 );
    end

    fprintf( '%d/%d     %f %f %f %f %f   %f %f %f %f    %d %d %d\n',...
        itrial, Ntrial, mysepar, mysepareeg, myseparnirschan, myseparnirssrc, myseparnirsdet, ...
        distmaxeeg, distmaxnirschan,distmaxnirssrc, distmaxnirsdet,...
        ns, nd, ne);
    itrial = itrial + 1;

end


%% output plot
flag_print = 0;

for ii = 1:length(rrrlist);
    rrr = rrrlist{ii};
    cc =        cclist{ii};
    ccsrc   =    ccsrclist{ii} ;
    ccdet    =   ccdetlist{ii}  ;
    mysepar  =     myseparlist(ii)  ;
    myseparnirschan   =    myseparnirschanlist(ii)  ;
    myseparnirssrc   =   myseparnirssrclist(ii) ;
    myseparnirsdet =      myseparnirsdetlist(ii)  ;

    fh = figure(330); clf; set(fh, 'color', 'white');  hold on;
    % plot nirs chans
    for kk = 1: length(cc),
        pnirschan = plot( cc(kk, 1), cc(kk, 2), 'markersize',10,'marker','s','markerfacecolor', [nirschangry],'markeredgecolor',[nirschangry]);
        plot( [ccsrc(kk,1) ccdet(kk,1)],  [ccsrc(kk,2) ccdet(kk,2)] ,'color',  [nirschangry] );
    end
    % plot eeg
    for kk = 1: numpos,
        if 1==rrr(kk), colr = srcolor;
        elseif 2==rrr(kk), colr = detcolor;
        elseif 3==rrr(kk), colr = eegcolor;
        else, error('should not be here'); end;
        if 1==rrr(kk),
            psrc = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
        elseif 2==rrr(kk),
            pdet = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
        elseif 3==rrr(kk), % eeg
            peeg = plot( xx(kk, 1), xx(kk, 2), 'markersize',8,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
        end
        text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', gry, 'fontsize',8 );
    end
    axis([-.2 6 -.2 5])
    box off
    axis off
    legend( [pnirschan psrc pdet peeg], 'NIRS Channel', 'Src', 'Det','EEG', 'fontsize',9 );
    legend boxoff;
    
    titlestr = [ ...
        num2str( (ii)  ) ' of ' ...
        num2str( length(rrrlist) ) ':  \Delta=' ...
        num2str( myseparlist(ii)  ) '   \Delta_{E}=' ...
        num2str( mysepareeglist(ii)  ) '   \Delta_{NC}=' ...
        num2str( myseparnirschanlist(ii) )  '   \Delta_{NS}=' ...
        num2str( myseparnirssrclist(ii) )  '   \Delta_{ND}=' ...
        num2str( myseparnirsdetlist(ii) ) ];



    title( titlestr, 'fontsize',9 );
    
    fstrn = [ 'fig_montage_'  num2str( ii ) '_'];
    set(gcf,'PaperPositionMode','auto');
    if flag_print, saveas(gcf, [ fstrn  '.png'  ] ); end
    pause(.2);

end













