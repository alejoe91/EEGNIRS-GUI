clear
figspath = 'C:\Users\ahmet\Documents\MATLAB\matlab-work\svn_matlab_projs\NIRSEEG\figs\';
fontsz = 22; fontsztic = 16; linsz = 3;

xminlim = 0;
xmaxlim = 180;
yminlim = 0;
ymaxlim = 1.04;
%%
load ccchistory_eeg_vft.mat
tccchistory_vft =  tccchistory;
ccchistory_vft = ccchistory;

load ccchistory_eeg_eyesclose1.mat
tccchistory_rs =  tccchistory;
ccchistory_rs = ccchistory;

figure(180);clf;
plot( tccchistory_vft, ccchistory_vft,'linewidth',linsz,'color','b'); hold on;
plot( tccchistory_rs, ccchistory_rs,'linewidth',linsz,'color','k','linestyle','--'); hold on;

axis(  [xminlim  xmaxlim  yminlim  ymaxlim]  );

xlabel('Time (s)', 'fontsize', fontsz); ylabel('Correlation', 'fontsize', fontsz);
set( gca, 'fontsize', fontsztic );

set(gca,'xTick',[xminlim:50:xmaxlim]);
set(gca,'yTick',[yminlim:.2:ymaxlim]);

set( gcf, 'PaperPositionMode', 'auto')
print(  '-dtiff', '-zbuffer', '-r600',  [figspath  'fig_corrhistory_eeg'] );
%%
load ccchistory_nirs_vft.mat
tccchistory_vft =  tccchistory;
ccchistory_vft = ccchistory;

load ccchistory_nirs_eyesclose1.mat
tccchistory_rs =  tccchistory;
ccchistory_rs = ccchistory;

figure(280);clf;
plot( tccchistory_vft, ccchistory_vft,'linewidth',linsz,'color','b'); hold on;
plot( tccchistory_rs, ccchistory_rs,'linewidth',linsz,'color','k','linestyle','--'); hold on;

axis(  [xminlim  xmaxlim  yminlim  ymaxlim]  );

xlabel('Time (s)', 'fontsize', fontsz);
%ylabel('Correlation', 'fontsize', fontsz);
set( gca, 'fontsize', fontsztic );

set(gca,'xTick',[xminlim:50:xmaxlim]);
set(gca,'yTick',[yminlim:.2:ymaxlim]);
set( gcf, 'PaperPositionMode', 'auto')
print(  '-dtiff', '-zbuffer', '-r600',  [figspath  'fig_corrhistory_nirs_oyx'] );
%%
load ccchistory_nirsR_vft.mat
tccchistory_vft =  tccchistory;
ccchistory_vft = ccchistory;

load ccchistory_nirsR_eyesclose1.mat
tccchistory_rs =  tccchistory;
ccchistory_rs = ccchistory;

figure(380);clf;
plot( tccchistory_vft, ccchistory_vft,'linewidth',linsz,'color','b'); hold on;
plot( tccchistory_rs, ccchistory_rs,'linewidth',linsz,'color','k','linestyle','--'); hold on;

axis(  [xminlim  xmaxlim  yminlim  ymaxlim]  );

xlabel('Time (s)', 'fontsize', fontsz);
%ylabel('Correlation', 'fontsize', fontsz);
set( gca, 'fontsize', fontsztic );

set(gca,'xTick',[xminlim:50:xmaxlim]);
set(gca,'yTick',[yminlim:.2:ymaxlim]);

set( gcf, 'PaperPositionMode', 'auto')
print(  '-dtiff', '-zbuffer', '-r600',  [figspath  'fig_corrhistory_nirs_deoxy'] );






