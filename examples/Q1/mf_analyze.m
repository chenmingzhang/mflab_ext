%% Analyzing output of the model

%% -----------------------------initialization-----------------------------
clear all;close all;
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat



%% ----------------------read calculation results -------------------------
% H is a struct to store hydraulic head calculated by MODFLOW.
% the number of outputs is determined by NTSP and INCODE in PER tab of excel
%   file.
% For example, if NTSP=100 and INCODE=10 are prescribed for three stree 
%   periods, then there will be 30 outputs (e.g., size(H)=30,1)
%   H(10).values are the hydraulic head at the end of the first stress period
%   H(21).values are the hydraulic head at the beginning of the third stress
%   period.
% H(i).totim is the real time [T] of the simulation

H=readDat([basename,'','.HDS']); % read the unformatted head file
B=readBud([basename, '.BGT']); % read the Budget file

% B is a struct storing the volumetric flow rate (L3/T) across each cell

%% ---------------------- a glance at the calculation domain---------------
gr.plotMesh('faceAlpha',0.15);view(3);


%% ---------------------plot hydraulic head information -------------------
a.fig=figure;
a.fs=12;
a.lw=1;
a.sub1=subplot(2,1,1);
% Pay attention to the data structure H here
% and also understand why int2str is requred
i=1 ;plot(gr.xm,H(i).values(:,:,1),'r-o','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on 
i=5 ;plot(gr.xm,H(i).values(:,:,1),'g-v','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on
i=10;plot(gr.xm,H(i).values(:,:,1),'b-*','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on
i=15;plot(gr.xm,H(i).values(:,:,1),'c-s','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on
i=20;plot(gr.xm,H(i).values(:,:,1),'m-x','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on
i=25;plot(gr.xm,H(i).values(:,:,1),'y-+','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on
i=30;plot(gr.xm,H(i).values(:,:,1),'k-d','Displayname',[int2str(H(i).totim),' min'],'LineWidth',a.lw);hold on
xlabel('x (m)','FontSize',a.fs,'FontWeight','bold')
ylabel('z (m)','FontSize',a.fs,'FontWeight','bold')
axis([0 100 -8 0] )
title('Hydraulic head profile');
legend('show','Location','Southeast')


% plot hydraulic head over time at a fixed point
a.sub2=subplot(2,1,2);
j=2;plot([H.totim],arrayfun(@(y) y.values(1,j,1),H),'r-o','Displayname',[int2str(gr.xm(j)),' m'],'LineWidth',a.lw) ;hold on
j=5;plot([H.totim],arrayfun(@(y) y.values(1,j,1),H),'g-v','Displayname',[int2str(gr.xm(j)),' m'],'LineWidth',a.lw) ;hold on
j=8;plot([H.totim],arrayfun(@(y) y.values(1,j,1),H),'b-*','Displayname',[int2str(gr.xm(j)),' m'],'LineWidth',a.lw) ;hold on
xlabel('Time (mim)','FontSize',a.fs,'FontWeight','bold')
ylabel('Hydraulic Head (m)','FontSize',a.fs,'FontWeight','bold')
title('Hydraulic head profile');
axis([0 1800 -8 0] )
legend('show','Location','East')


print(a.fig,'head.png','-dpng')   % png figure output
print(a.fig,'head.tif','-dtiff','-r70')   % tif figure output
print(a.fig,'head.eps','-depsc')  % vectorized figure output (publishing standard)
                               % it can be viewed using ghostscript, evince
			       % or by insterting it into word
saveas(a.fig,'head.fig','fig') % save the matlab figure

%%-----------------------------plot head contour -----------------------
b.fig=figure;
b.fs=12;
b.lw=1;
b.sub1=subplot(4,1,1);
i=1 ;contourf(squeeze(gr.XM)' , squeeze(gr.ZM)' , squeeze(H(i).values)')
caxis([-6,0]);colorbar
b.sub2=subplot(4,1,2);
i=10;contourf(squeeze(gr.XM)' , squeeze(gr.ZM)' , squeeze(H(i).values)')
caxis([-6,0]);colorbar
b.sub3=subplot(4,1,3);
i=20;contourf(squeeze(gr.XM)' , squeeze(gr.ZM)' , squeeze(H(i).values)')
caxis([-6,0]);colorbar
b.sub4=subplot(4,1,4);
i=30;contourf(squeeze(gr.XM)' , squeeze(gr.ZM)' , squeeze(H(i).values)')
caxis([-6,0]);colorbar
print(b.fig,'coutour.png','-dpng')   % png figure output
print(b.fig,'coutour.tif','-dtiff','-r70')   % tif figure output
print(b.fig,'coutour.eps','-depsc')  % vectorized figure output (publishing standard)
                               % it can be viewed using ghostscript, evince
			       % or by insterting it into word
saveas(b.fig,'countour.fig','fig')  % save the matlab figure
