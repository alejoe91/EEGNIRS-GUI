clear;

N =4;
Nt = 100;
Fs = 250;
tt = [1:Nt]./Fs;

xxe = randn(Nt,N).*.2;
for kk=1:N,
    xxe(:,kk) = xxe(:,kk) + sin(2.*pi.*kk.*tt)';
end

figure(100)
for kk=1:N,
    subplot(N,1,kk), plot(tt, xxe(:,kk))
end

%cove = cov(xxe);
cove = xxe' * xxe;

[EE,EV] = eig(cove);    lame = diag(EV);     [lame,indx] = sort(lame,1,'descend');    evece = EE(:,indx);

aae = xxe * evece;
figure(200)
for kk=1:N,
    std(aae(:,kk))
    subplot(N,1,kk), plot(tt, aae(:,kk)); axis([tt(1) tt(end) -2 2])
end

for ii=1:N,
    for jj=ii:N,
        fprintf('%d %d %f\n',ii,jj,         sum( aae(:,ii).*aae(:,jj) ) );
    end
end

