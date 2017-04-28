clear; fprintf('\n');

%% parameters

% biophysical params
delta = .04;% diffusive step size
pa = .01;% probabiliy of absorption at each time step

% geometry params
xdet = 1.; % detector position
radiusdet = .05;% detector radius

% comput params
maxnumphoton = 100000;
maxiter = 100000000;

%% calcs

countdetected = 0;
for iphoton = 1:maxnumphoton,
    
    nwrite = ceil(maxnumphoton./10);
    if 0==mod(iphoton,nwrite), fprintf('photon %d / %d, detected %d\n',iphoton,maxnumphoton,countdetected);  end

    % initialize position
    x = zeros(1,2);
    isabsorbed = false; isescaped = false; isdetected = false;
    xpath = [];
    countiter = 1;
    
    % iterate time steps updating status of photon
    while countiter<=maxiter && ~isabsorbed && ~isescaped && ~isdetected,
                
        xpath(countiter,:) = x;
        countiter = countiter + 1;
        
        % scatter
        x = x + delta.*randn(1,2);
        % absorb
        if rand(1,1)<pa, isabsorbed = true; end;
        % escape domain
        if 0<x(2), isescaped = true; end;
        % detected
        if xdet-radiusdet<=x(1) & x(1)<=xdet+radiusdet & 0<x(2), isdetected = true; end

    end

    if isdetected,
        xpath(countiter,:) = x;
        countdetected = countdetected + 1;
        xpaths{countdetected} = xpath;
       
        if 0
            % plot path of one photon
            figure(320);clf; hold on; grid on;
            plot(xpath(:,1),xpath(:,2),'g-');
            plot(xpath(1,1),xpath(1,2),'r.');
            plot(xpath(end,1),xpath(end,2),'b.');
            title([ 'Photon# ' num2str(iphoton) ] );
            drawnow
            pause(1)
        end

    end

end
fprintf('Detected: %d / %d\n',countdetected,maxnumphoton);

%% plot paths
if 0<countdetected,

    figure(210);clf; hold on; grid on;
    for ii = 1:length(xpaths),
        xpath = xpaths{ii};
        plot(xpath(:,1),xpath(:,2),'g-');
        plot(xpath(1,1),xpath(1,2),'r.');
        plot(xpath(end,1),xpath(end,2),'b.');
    end
    title([ num2str(countdetected) '/' num2str(maxnumphoton) ] );
    xlabel('x'); ylabel('y');
    
    % calc mean
    xpath = [];
    for ii = 1:length(xpaths),
        xpath = [xpath; xpaths{ii}];
    end
    xmean = mean(xpath(:,1));
    ymean = mean(xpath(:,2));
    fprintf( 'COG: %f %f\n',xmean,ymean);
    plot( xmean,ymean,'k*','markersize',12); plot( xmean,ymean,'ko','markersize',12);
   
end