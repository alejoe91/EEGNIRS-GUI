%% PCA: for X where each row is an observation, and each column a variable,
% COV(X) is the covariance matrix.
if 0
    Q = cov( xxfn ); [EE,EV] = eig(Q); ev = diag(EV); [ev,indx] = sort(ev,1,'descend');
    EE = EE(:,indx); figure(220); clf; plot( ev,'.-' );
    aa = xxfn * EE;
    %% VIEW PCA MODES
    nbeg = Fs.*19; nend = nbeg + Fs.*20;
    tt = [nbeg:nend]./Fs;
    imode = 5;
    figure(420); clf; plot( tt, aa(nbeg:nend, imode) );
end
%% ICA
if 0,
    % Row j of W is the jth eigenvector.
    % s=Wx and x=As
    [icasig, A, W] = fastica ( xxfn' );
    icasig = icasig';
    %% VIEW ICA
    nbeg = Fs.*19; nend = nbeg + Fs.*20;
    tt = [nbeg:nend]./Fs;
    figure(550); clf;
    for imode = 1:numchan,
        plot( tt, icasig(nbeg:nend, imode) ); title( num2str(imode) );
        pause ( 1 );
    end
    %%
    iaa = xxfn * W';
    %%
    figure(560); clf;
    imode = 2,
    plot( tt, iaa(nbeg:nend, imode) ); title( num2str(imode) );
    %% REMOVE EYEBLINK
    imoderemove = [2, 10];
    Wc = W;
    Wc(imoderemove, :) = 0;
    xxfnc = ( xxfn * Wc' ) * A';
    %% VIEW CLEAN SIGNAL
    ichan = 2;
    figure(660); clf; hold on;
    plot( tt, xxfn(nbeg:nend, ichan), 'b' );
    plot( tt, xxfnc(nbeg:nend, ichan), 'r' );
end