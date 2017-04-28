clear;
%{
Distributes
sensors over whole head
in the 130 chan arrangement, there are 10 in midline, 60 on one side
%}
%% params
irandseed = 55;
rand('seed', irandseed );
writefreq = 100;
flag_print = 0;
maxsepar = 1.5; % cutoff (largest allowed) src-det distance for nirs channel
Ntrial = 1e9;
Maxtrial = 1e12;
labeloffset = .15;
%% max numbers of sensors
ns = 24; % source 1
nd = 24; % detector 2
ne = 19; % eeg 3
fprintf('ne %d, ns %d, nd %d\n',ne, nd, ns);
%%
gry = [.6 .6 .6];
eeggry = [.4 .0 .4];
nirschangry = [.6 .6 .6];
nirschanconnectgry  = [.9 .9 .9];
portlabelgry  = [.7 .7 .7];
srcolor = 'y';
detcolor = 'g';
eegcolor = eeggry;
%% input labels and their positions
fid = fopen('positions_full.txt');
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

if 0
%% insert additional positions on right side
[xx, label] = insertadditionalpositions( xx, label );

%% insert left side by taking the negative of x axis
[xx, label] = insertleftpositions( xx, label );
end

%% 2d plot
if 1,
    figure(1010);clf; hold on;
    labeloffset = .2;
    for kk = 1:length(label),
        plot( xx(kk, 1), xx(kk, 2), 'ro' );
        text( xx(kk, 1), xx(kk, 2)+labeloffset, label( kk ), 'color', gry, 'fontsize',7 );
    end
    plot( [-5.5 5.5],[0 0], 'color',gry,'linestyle',':'); plot( [0 0],[-5.5 4.5], 'color',gry,'linestyle',':')
%    axis([-6 6 -6 5])
end
%% 3d plot
if 1,
    figure(1020);clf; hold on;
    for kk = 1:length(label),
        if xx(kk,1)<0, continue; end;
        plot3( xx(kk, 1), xx(kk, 2), xx(kk, 3), 'ro' );
        text( xx(kk, 1)+labeloffset, xx(kk, 2), xx(kk, 3), label( kk ), 'color', gry );
    end
    %    axis([-1 14 -1 14 -1 14])   
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

    rrr = zeros( 1, numpos );

    % choose eeg sensor positions
    posindxs = [1:numpos];
    pop = posindxs;
    eegpos = randsample( pop, ne );
    rrr( eegpos ) = 3;
    pop = setdiff( pop, eegpos); % remove eeg positions from candidates

    % find indx of bottom half
    ss = [];
    for kk = 1: numpos,
        if 0>xx(kk,2), ss = [ss; kk]; end
    end
    pop = setdiff( pop, ss); % remove bottom half from candidates

    % choose nirs source positions
    nsmax = min(ns, length(pop) );
    nirsspos = randsample( pop, nsmax );
    rrr( nirsspos ) = 1;
    pop = setdiff( pop, nirsspos); % remove nirs src positions

    % choose nirs detector positions
    ndmax = min(nd, length(pop) );
    nirsdpos = randsample( pop, ndmax );
    rrr( nirsdpos ) = 2;


    %% set nirs channels
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
            elseif 0==rrr(kk), colr = 'r';
            else, error('should not be here'); end;
            if 1==rrr(kk), % nirs src
                psrc = plot( xx(kk, 1), xx(kk, 2), 'markersize',6,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 2==rrr(kk),
                pdet = plot( xx(kk, 1), xx(kk, 2), 'markersize',6,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 3==rrr(kk), % eeg
                peeg = plot( xx(kk, 1), xx(kk, 2), 'markersize',10,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
            elseif 0==rrr(kk),
                plot( xx(kk, 1), xx(kk, 2), 'markersize',3,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
            else, error('should not be here'); end;
            text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', gry, 'fontsize',8 );
        end
        axis([-6 6 -6 7])
        box off
        axis off
        if numnirschan, legend( [pnirschan psrc pdet peeg], 'NIRS Channel', 'Src', 'Det','EEG','location','southwest'  );
        else, legend( [psrc pdet peeg],  'Src', 'Det','EEG' ,'location','southwest' ); end;
    end
    %% extract eeg sensors
    indx = find( 3==rrr );
    xxeeg = xx( indx, : );

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
    iii = 1;
    for kk = 1:numpos,
        if 0<xx(kk,2), % must include only frontal positions
            distmin = 10000;
            for nn = 1:length( cc ),
                dist = sqrt( (xx(kk,1)-cc(nn,1)).^2 + (xx(kk,2)-cc(nn,2) ).^2 );
                if dist < distmin,
                    distmin = dist;
                end
            end
            distnirs(iii) = distmin;
            iii = iii + 1;
        end
    end
    % largest distance (from any position) to the nearest nirs chan
    distmaxnirschan = max( distnirs );

    %% determine distance to nearest nirs source
    iii = 1;
    for kk = 1:numpos,
        if 0<xx(kk,2), % must include only frontal positions
            distmin = 10000;
            for nn = 1:numpos,
                if 1==rrr(nn) & kk~=nn,
                    dist = sqrt( (xx(kk,1)-xx(nn,1)).^2 + (xx(kk,2)-xx(nn,2) ).^2 );
                    if dist < distmin,
                        distmin = dist;
                    end
                end
            end
            distnirs(iii) = distmin;
            iii = iii + 1;
        end
    end
    % largest distance (from any position) to the nearest nirs chan
    distmaxnirssrc = max( distnirs );
    %% determine distance to nearest nirs det
    iii = 1;
    for kk = 1:numpos,
        if 0<xx(kk,2), % must include only frontal positions
            distmin = 10000;
            for nn = 1:numpos,
                if 2==rrr(nn) & kk~=nn,
                    dist = sqrt( (xx(kk,1)-xx(nn,1)).^2 + (xx(kk,2)-xx(nn,2) ).^2 );
                    if dist < distmin,
                        distmin = dist;
                    end
                end
            end
            distnirs(iii) = distmin;
            iii = iii + 1;
        end
    end
    % largest distance (from any position) to the nearest nirs chan
    distmaxnirsdet = max( distnirs );
    %%

    % crietrion
    myvec = [distmaxeeg distmaxnirschan distmaxnirssrc distmaxnirsdet];
    ddd = sqrt( mean(myvec.^2) );
    % update minima
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
        myne(countgood) = ne;
        mynumnirschan(countgood) = numnirschan;
        myndmax(countgood) = ndmax;
        mynsmax(countgood) = nsmax;

        %% good plot
        fh = figure(220); clf; set(fh, 'color', 'white');  hold on;
        %% plot nirs chans
        for kk = 1: length(cc), % connectors
            plot( [ccsrc(kk,1) ccdet(kk,1)],  [ccsrc(kk,2) ccdet(kk,2)] ,'color',  [nirschanconnectgry] );
        end
        for kk = 1: length(cc), % chans
            pnirschan = plot( cc(kk, 1), cc(kk, 2), 'markersize',10,'marker','s','markerfacecolor', [nirschangry],'markeredgecolor',[nirschangry]);
        end
        %% plot eeg
        for kk = 1: numpos,
            if 1==rrr(kk), colr = srcolor;
            elseif 2==rrr(kk), colr = detcolor;
            elseif 3==rrr(kk), colr = eegcolor;
            elseif 0==rrr(kk), colr = 'r';
            else, error('should not be here'); end;
            if 1==rrr(kk),
                psrc = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 2==rrr(kk),
                pdet = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
            elseif 3==rrr(kk), % eeg
                peeg = plot( xx(kk, 1), xx(kk, 2), 'markersize',8,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
            elseif 0==rrr(kk),
                plot( xx(kk, 1), xx(kk, 2), 'markersize',3,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
            end
            text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', gry, 'fontsize',8 );
        end
        axis([-4.8  5. -5 5.3])
        box off
        axis off
        legend( [pnirschan psrc pdet peeg], 'NIRS Channel', 'Src', 'Det','EEG','location','southwest' );
        legend boxoff;
        titlestr1 = [ ...
            num2str( countgood ) ' (' ...
            num2str( itrial  ) ' of ' ...
            num2str( Ntrial ) '):  \Delta=' ...
            num2str( mysepar ) '  \Delta_{E}=' ...
            num2str( mysepareeg  ) '   \Delta_{c}=' ...
            num2str( myseparnirschan )  '   \Delta_{s}=' ...
            num2str( myseparnirssrc )  '   \Delta_{d}=' ...
            num2str( myseparnirsdet ) ];
        titlestr2 =  [ ...
            'N_E=' num2str(ne)...
            ', N_c=' num2str(numnirschan)...
            ', N_s=' num2str(nsmax)...
            ', N_d=' num2str(ndmax)...
            ];
        title( {titlestr1; titlestr2}, 'fontsize',9 );

        if flag_print,
            fstrn = [ 'fig_montage_'  'randseed_' num2str(irandseed) '_' num2str( countgood ) '_'];
            set(gcf,'PaperPositionMode','auto');
            if flag_print, saveas(gcf, [ fstrn  '.png'  ] ); end
        end
        
        countgood = countgood + 1;
        pause(.2);
    end
    if 0==mod(itrial,writefreq),
        fprintf( '%d/%d     %d %d %d %d   %f    %f %f %f %f   %f %f %f %f \n',...
            itrial, Ntrial, ne,nsmax,ndmax,numnirschan,...
            mysepar, mysepareeg, myseparnirschan, myseparnirssrc, myseparnirsdet, ...
            distmaxeeg, distmaxnirschan,distmaxnirssrc, distmaxnirsdet...
            );
    end
    
    itrial = itrial + 1;
end
return;

%% output plot
for countgood = 1:length(rrrlist),
    rrr = rrrlist{countgood} ;
    cc = cclist{countgood} ;
    ccsrc =  ccsrclist{countgood} ;
    ccdet = ccdetlist{countgood} ;
    mysepar = myseparlist(countgood) ;
    mysepareeg = mysepareeglist(countgood) ;
    myseparnirschan = myseparnirschanlist(countgood) ;
    myseparnirssrc = myseparnirssrclist(countgood) ;
    myseparnirsdet = myseparnirsdetlist(countgood) ;
    ne = myne(countgood) ;
    numnirschan = mynumnirschan(countgood) ;
    ndmax = myndmax(countgood) ;
    nsmax = mynsmax(countgood) ;

    %  plot channels
    fh = figure(330); clf; set(fh, 'color', 'white');  hold on;
    % plot nirs chans
    for kk = 1: length(cc), % channel connectors
        plot( [ccsrc(kk,1) ccdet(kk,1)],  [ccsrc(kk,2) ccdet(kk,2)] ,'color',  [nirschanconnectgry] , 'linestyle','-');
    end
    for kk = 1: length(cc), % channels
        pnirschan = plot( cc(kk, 1), cc(kk, 2), 'markersize',10,'marker','s','markerfacecolor', [nirschangry],'markeredgecolor',[nirschangry]);
    end
    % plot  src, det, eeg
    for kk = 1: numpos,
        if 1==rrr(kk), colr = srcolor;
        elseif 2==rrr(kk), colr = detcolor;
        elseif 3==rrr(kk), colr = eegcolor;
        elseif 0==rrr(kk), colr = 'r';
        else, error('should not be here'); end;
        if 1==rrr(kk),
            psrc = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
        elseif 2==rrr(kk),
            pdet = plot( xx(kk, 1), xx(kk, 2), 'markersize',5,'marker','o','markerfacecolor', [colr],'markeredgecolor',[colr]);
        elseif 3==rrr(kk), % eeg
            peeg = plot( xx(kk, 1), xx(kk, 2), 'markersize',8,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
        elseif 0==rrr(kk),
            plot( xx(kk, 1), xx(kk, 2), 'markersize',3,'marker','o','markerfacecolor', [colr],'markeredgecolor','none');
        end
        text( xx(kk, 1)+labeloffset, xx(kk, 2), label( kk ), 'color', portlabelgry, 'fontsize',8 );
    end
    axis([-4.8  5. -5 5.3])
    box off
    axis off
    legend( [pnirschan psrc pdet peeg], 'NIRS Channel', 'Src', 'Det','EEG','location','southwest' );
    legend boxoff;
    titlestr1 = [ ...
        num2str( countgood  ) ' of ' ...
        num2str( length(rrrlist) ) ':  \Delta=' ...
        num2str( mysepar ) '  \Delta_{E}=' ...
        num2str( mysepareeg  ) '   \Delta_{c}=' ...
        num2str( myseparnirschan )  '   \Delta_{s}=' ...
        num2str( myseparnirssrc )  '   \Delta_{d}=' ...
        num2str( myseparnirsdet ) ];
    titlestr2 =  [ ...
        'N_E=' num2str(ne)...
        ', N_c=' num2str(numnirschan)...
        ', N_s=' num2str(nsmax)...
        ', N_d=' num2str(ndmax)...
        ];
    title( {titlestr1; titlestr2}, 'fontsize',9 );
    pause(1);
end




