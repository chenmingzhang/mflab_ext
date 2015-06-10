%% Analyzing output of the model
% TO 091011 091129 120413
 
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

movie1 = false;
movie2 = false;
movie3 = false;

Rkm  = R/1000;
R1km = R1/1000;

Hds=readDat(  [basename,'','.HDS']); % read the unformatted head file
hrange = ContourRange(Hds,50); % get contour intervals

B= readBud6([basename, '.BGT']); % read the Budget file
if gr.Nx>1
    B = mf_Psi(gr,B,1:gr.Ny);
    prange = ContourRange(B,50,'Psi');
else
    movie1=false;
end

%%

time   = [Hds.time];

[PERnams,PERvals] = getPeriods(basename,'PER',{'year','month','day'});
date   = datenum(PERvals);

timeYr = time./aYear; 

sunny = sun([B.period]); sunPos = [0.03 0.85 0.08 0.08];

% R1 comes from underneath

fsz = 14;

%% Setup movie1

if movie1

    figure('pos',screenPos(0.75));
    ax1 = subplot(2,2,1,'nextplot','add','ylim',[-5          2],'xlim',[0 R1km],'fontSize',fsz);
    ax2 = subplot(2,2,3,'nextplot','add','ylim',[gr.zGr(end),0],'xlim',[0 R1],'fontSize',fsz);
    ax3 = subplot(2,2,2,'nextplot','add','ylim',[-5          2],'xlim',[0 R1km],'fontSize',fsz);
    ax4 = subplot(2,2,4,'nextplot','add','ylim',[gr.zGr(end),0],'xlim',[0 R1],'fontSize',fsz);

    %% Use zone budget to get budget overview

    % zonebudget(B(end));

    %axWeather = axes('pos',[0.49,0.45,0.1,0.1],'nextPlot','add','color','none','visible','off','yDir','reverse');
    axWeather = axes('pos',sunPos,'nextPlot','add','color','none','visible','off','yDir','reverse');
    Sun  = imread('Sun.png');
    Rain = imread('Rain','png');
    hSun = image(Sun,'parent',axWeather);
    hRain= image(Rain,'parent',axWeather);

    %% show dynamically varying water levels and stream lines

    vidObj = VideoWriter([basename sprintf('_scen_%d',scenario)]);
    vidObj.FrameRate = 15;
    vidObj.Quality = 80;
    vidObj.open();

    leg = {   'layer 1',   'layer 2',   'layer 3'};

    for it=1:numel(timeYr)
        s = sprintf('Head and stream lines Terwisscha, scenario %d, t=%.2f y',scenario,timeYr(it));
        if it==1
            ht1  = title(ax1,s);
            ht2  = title(ax3,s);

            hp1 = plot(ax1,gr.xkm,Hds(it).values(1,:,1),'r');
            hp2 = plot(ax1,gr.xkm,Hds(it).values(1,:,2),'b');
            if gr.Nlay>2
                hp3 = plot(ax1,gr.xkm,Hds(it).values(1,:,end),'k');
            end

            hp4 = plot(ax3,gr.xkm,Hds(it).values(2,:,1),'r');
            hp5 = plot(ax3,gr.xkm,Hds(it).values(2,:,2),'b');
            if gr.Nlay>2
                hp6 = plot(ax3,gr.xkm,Hds(it).values(2,:,end),'k');
            end

            plot(ax1,[R,gr.xGr(end)]/1000,[0 0],'g','linewidth',2);
            text(Rkm,0,'Agricultural land -->','VerticalAlignment','Bottom','fontSize',12,'parent',ax1);

            legend(ax1,leg{1:gr.Nlay},4);

            plot(ax3,[R,gr.xGr(end)]/1000,[0 0],'g','linewidth',2);
            text(Rkm,0,'Agricultural land -->','VerticalAlignment','Bottom','fontSize',12,'parent',ax3);

            legend(ax3,leg{1:gr.Nlay},4);

            gr.plotXSec(ax2,1,'all','hlines','k','edgecolor',grey,'faceColor','none','title','streamlines');
            gr.plotXSec(ax4,2,'all','hlines','k','edgecolor',grey,'faceColor','none','title','streamlines');

            xlabel(ax1,'r [km]');
            xlabel(ax3,'r [km]');

            ylabel(ax1,'head [m]');
            ylabel(ax3,'head [m]');
            ylabel(ax2,'elevation [m]');
            ylabel(ax4,'elevation [m]');

            [~,hp7] = contour(ax2,gr.xp,gr.zp,B(it).Psi(:,:,1),prange);
            [~,hp8] = contour(ax4,gr.xp,gr.zp,B(it).Psi(:,:,2),prange);

            plot(ax2,[R,gr.xGr(end)],[0 0],'g','linewidth',2);
            text(R,0,'Agricultural land -->','VerticalAlignment','Top','fontSize',12,'parent',ax2);  

            plot(ax4,[R,gr.xGr(end)],[0 0],'g','linewidth',2);
            text(R,0,'Agricultural land -->','VerticalAlignment','Top','fontSize',12,'parent',ax4);        

        else
            set(ht1,'string',s);
            set(ht2,'string',s);

            set(hp1,'ydata',Hds(it).values(1,:,1));
            set(hp2,'ydata',Hds(it).values(1,:,2));
            if gr.Nlay>2
                set(hp3,'ydata',Hds(it).values(1,:,end));
            end

            set(hp4,'ydata',Hds(it).values(2,:,1));
            set(hp5,'ydata',Hds(it).values(2,:,2));        
            if gr.Nlay>2
                set(hp6,'ydata',Hds(it).values(2,:,end));
            end

            set(hp7 ,'zdata',B(it).Psi(:,:,1));
            set(hp8 ,'zdata',B(it).Psi(:,:,2));

        end

        if sunny(it);
            set(hSun,'visible','on'); set(hRain,'visible','off');
        else
            set(hSun,'visible','off'); set(hRain,'visible','on');
        end

        vidObj.writeVideo(getframe(gcf));   
    end
    vidObj.close();

end

%% Dynamically show the drawdown, hat is the fluctuations combined, and the drawdown in two subplots

if movie2
    
    figure('pos',screenPos(0.75));
    ax5 = subplot(3,1,1,'nextplot','add','ylim',[-5   2],'xlim',[0 R1]/1000,'fontSize',fsz);
    ax6 = subplot(3,1,2,'nextplot','add','ylim',[1e-3 10],'ytick',[1e-3 1e-2 1e-1 1e0 10],'xlim',[0 R1]/1000,'yscale','log','yDir','reverse',...
        'yGrid','on','yMinorGrid','on','fontSize',fsz);
    ax7 = subplot(3,1,3,'nextplot','add','ylim',[0    3],'xlim',[0 R1]/1000,'yscale','lin','yDir','reverse',...
        'yGrid','on','yMinorGrid','on','fontSize',fsz);

    axWeather = axes('pos',sunPos,'nextPlot','add','color','none','visible','off','yDir','reverse');
    Sun  = imread('Sun.png');
    Rain = imread('Rain','png');
    hSun = image(Sun,'parent',axWeather);
    hRain= image(Rain,'parent',axWeather);

    vidObj = VideoWriter([basename '_dd' sprintf('_scen_%d',scenario)]);
    vidObj.FrameRate = 15;
    vidObj.Quality   = 80;
    vidObj.open();
    for it=1:numel(timeYr)
        s = sprintf('Drawdown Terwisscha, scenario %d, XSection, t=%.2f y',scenario,timeYr(it));
        if it==1
             ht1  = title(ax5,s);

             if gr.Nlay>2
                 hp3 = plot(ax5,gr.xkm,Hds(it).values(1,:,end),'k');
             end
             hp2 = plot(ax5,gr.xkm,Hds(it).values(1,:,2),'b');
             hp1 = plot(ax5,gr.xkm,Hds(it).values(1,:,1),'r');

             if gr.Nlay>2
                 hp6 = plot(ax5,gr.xkm,Hds(it).values(2,:,end),'k');
             end
             hp5 = plot(ax5,gr.xkm,Hds(it).values(2,:,2),'b');
             hp4 = plot(ax5,gr.xkm,Hds(it).values(2,:,1),'r');

             if gr.Nlay>2
                 hp9 = plot(ax6,gr.xkm,Hds(it).values(2,:,end)-Hds(it).values(1,:,end),'k');
             end
             hp8 = plot(ax6,gr.xkm,Hds(it).values(2,:,2)-Hds(it).values(1,:,2),'b');
             hp7 = plot(ax6,gr.xkm,Hds(it).values(2,:,1)-Hds(it).values(1,:,1),'r');

             if gr.Nlay>2
                 hp12 = plot(ax7,gr.xkm,Hds(it).values(2,:,end)-Hds(it).values(1,:,end),'k');
             end
             hp11 = plot(ax7,gr.xkm,Hds(it).values(2,:,2)-Hds(it).values(1,:,2),'b');
             hp10 = plot(ax7,gr.xkm,Hds(it).values(2,:,1)-Hds(it).values(1,:,1),'r');

    %          if scenario==1
    %              plot(ax5,gr.xkm, sDeGlee,'g','lineWidth',3); % De Glee
    %              plot(ax6,gr.xkm,-sDeGlee,'g','lineWidth',3); % De Glee
    %              plot(ax7,gr.xkm,-sDeGlee,'g','lineWidth',3); % De Glee
    %          end

             plot(ax5,[R,gr.xGr(end)]/1000,[0 0],'g','linewidth',2);
             text(Rkm,0,'Agricultural land -->','VerticalAlignment','Bottom','fontSize',12,'parent',ax5);

             plot(ax6,[R,gr.xGr(end)]/1000,[0.001 0.001],'g','linewidth',2);
             text(Rkm,0.001,'Agricultural land -->','VerticalAlignment','Top','fontSize',12,'parent',ax6);

             plot(ax7,[R,gr.xGr(end)]/1000,[0 0],'g','linewidth',2);
             text(Rkm,0,'Agricultural land -->','VerticalAlignment','Top','fontSize',12,'parent',ax7);

             legend(ax5,leg{1:gr.Nlay},4);
             legend(ax6,leg{1:gr.Nlay},4);
             legend(ax7,leg{1:gr.Nlay},4);

             xlabel(ax5,'r [km]');
             xlabel(ax6,'r [km]');
             xlabel(ax7,'r [km]');

             ylabel(ax5,'head [m]');
             ylabel(ax6,'drawdown [m]');
             ylabel(ax7,'drawdown [m]');

             plot(ax6,get(ax6,'xlim'),[0.01 0.01],'k','lineWidth',3);
             plot(ax6,get(ax6,'xlim'),[0.1  0.1 ],'k','lineWidth',3);
             text(7,0.01,'1 cm','parent',ax6,...
                 'fontsize',16,'fontWeight','bold','fontAngle','italic','verticalAlignment','bottom');
             text(7,0.1,'10 cm','parent',ax6,...
                 'fontsize',16,'fontWeight','bold','fontAngle','italic','verticalAlignment','bottom');

        else
            set(ht1,'string',s);
            set(hp1,'ydata',Hds(it).values(1,:,1));
            set(hp2,'ydata',Hds(it).values(1,:,2));
            set(hp4,'ydata',Hds(it).values(2,:,1));
            set(hp5,'ydata',Hds(it).values(2,:,2));
            set(hp7,'ydata',Hds(it).values(2,:,1)-Hds(it).values(1,:,1));        
            set(hp8,'ydata',Hds(it).values(2,:,2)-Hds(it).values(1,:,2)); 
            set(hp10 ,'ydata',Hds(it).values(2,:,1)-Hds(it).values(1,:,1));        
            set(hp11,'ydata',Hds(it).values(2,:,2)-Hds(it).values(1,:,2)); 
            if gr.Nlay>2
                set(hp3,'ydata',Hds(it).values(1,:,end));
                set(hp6,'ydata',Hds(it).values(2,:,end));
                set(hp9 ,'ydata',Hds(it).values(2,:,end)-Hds(it).values(1,:,end)); 
                set(hp12,'ydata',Hds(it).values(2,:,end)-Hds(it).values(1,:,end)); 
            end
        end

        if sunny(it);
            set(hSun,'visible','on'); set(hRain,'visible','off');
        else
            set(hSun,'visible','off'); set(hRain,'visible','on');
        end

        vidObj.writeVideo(getframe(gcf));   
    end
    vidObj.close();

end 

%% Show drawdowns as function of time

% Compute drawdown by subtracting the head of the parallelly computed axial 
% models
% Dd.values(1,:,:) is the drawdown of the DRN mdoel
% Dd.values(2,:,:) is the drawdown of the GHB model
Dd = Hds;  % drawdown
DB = B;  % water budget differences
Qdiff = NaN(numel(B),numel(B(1).label),3);

for it=1:numel(Hds)
    Dd(it).time  = Dd(it).time/365.24; % change to years
    Dd(it).values = NaN(gr.size(1)/2,gr.Nx,gr.Nlay);
    Dd(it).values(1,:,:) = Hds(it).values(2,:,:)-Hds(it).values(1,:,:); % DRN model Q = 7.5 Mm3/a
    Dd(it).values(2,:,:) = Hds(it).values(4,:,:)-Hds(it).values(3,:,:); % GHB model Q = 7.5 Mm3/a
    Dd(it).values(3,:,:) = Hds(it).values(6,:,:)-Hds(it).values(5,:,:); % DRN model Q =2.5 Mm3/a
    Dd(it).NROW=size(Dd(it).values,1);
    Dd(it).rows = 1:Dd(it).NROW;
    
    for j=1:numel(B(it).label)
        DB(it).term{j} = NaN(gr.size(1)/2,gr.Nx,gr.Nlay);
        DB(it).term{j}(1,:,:) = B(it).term{j}(1,:,:)-B(it).term{j}(2,:,:);
        DB(it).term{j}(2,:,:) = B(it).term{j}(3,:,:)-B(it).term{j}(4,:,:);
        DB(it).term{j}(3,:,:) = B(it).term{j}(5,:,:)-B(it).term{j}(6,:,:);
        Qdiff(it,j,:) = sum(sum(DB(it).term{j},2),3);
    end
end


% Put all drawdowns in a single array (Nlay,Nx,Nt)
drawd1 = NaN(gr.Nz,gr.Nx,numel(Dd));
drawd2 = NaN(gr.Nz,gr.Nx,numel(Dd));

for it = numel(Dd):-1:1
    drawd1(:,:,it) = XS(Dd(it).values(1,:,:));
    drawd2(:,:,it) = XS(Dd(it).values(2,:,:));
end

%% extract the highest and the lowest drawdowns in any one year
year = PERvals(:,strmatchi('year',PERnams));
yr   = unique(year);
yr = yr(2:end);

for iYr = numel(yr):-1:1
    drawd1Myr(:,:,iYr) = max(drawd1(:,:,year==yr(iYr)),[],3);
    drawd1myr(:,:,iYr) = min(drawd1(:,:,year==yr(iYr)),[],3);
    drawd2Myr(:,:,iYr) = max(drawd2(:,:,year==yr(iYr)),[],3);
    drawd2myr(:,:,iYr) = min(drawd2(:,:,year==yr(iYr)),[],3);
end
drawd1Myr = mean(drawd1Myr,3);
drawd1myr = mean(drawd1myr,3);
drawd2Myr = mean(drawd2Myr,3);
drawd2myr = mean(drawd2myr,3);

%% Plot drawdown as function of xm at end of simulation
figure; hold on; set(gca,'yDir','reverse');
plot(gr.xm,drawd1(:,:,end),'r');
plot(gr.xm,drawd2(:,:,end),'b');
legend('M1','M2',4);

%% plot drawdown as function of time for different distances from well

figure('name','dd vs time','pos',screenPos(1,0.75)); 
axes('nextPlot','add','yDir','reverse','yGrid','off','fontSize',fsz);

% desired distances
x = [2500 3500 4500 5500 6500 7500];
% replace x by closest distances and get the x-index of those distances
% so we can pick the x-coordinate from gr.xm(ix)
for ix=numel(x):-1:1
    Ix(ix) = find(gr.xm<=x(ix),1,'last');
    x(ix)=gr.xm(Ix(ix));
end

% plot results, for only plot one of the four graphs, just to get the legend
k=1; leg={};
for ix=Ix
    plot(date,squeeze(drawd1(1,ix,:)),[mf_color(k),'--'],'lineWidth',2);
    leg{k} = sprintf('x=%.0fm',gr.xm(ix)); %#ok
    k=k+1;
end

% then plot the other curvecs for the distances using the same color as the
% first
k=1;
for ix=Ix
    plot(date,squeeze(drawd1(2,ix,:)),[mf_color(k),'-' ],'lineWidth',2);
    plot(date,squeeze(drawd2(1,ix,:)),[mf_color(k),'--'],'lineWidth',1);
    plot(date,squeeze(drawd2(2,ix,:)),[mf_color(k),'-' ],'lineWidth',1);
    k=k+1;
end

legend(leg{:},3); datetick(); xlabel('time'); ylabel('drawdown [m]');
title('Drawdown versus time at various distances from pumping station');

%% Checking the computed heads at the locations of the piezoms1

if 0  % takes to much time, skip
    figure('pos',screenPos(0.75));
    ax = axes('nextPlot','add','yDir','normal','fontSize',fsz,'xGrid','on','yGrid','on','ylim',[0 1.2]);
    xlabel(ax,'time [d] -->'); ylabel('Heads [m]');
    title(sprintf('Heads at differnt R vs time, scenario %d',scenario))
    xlabel('time [y] -->');
    ylabel('Head [m]');

    leg={};
    for ip=1:numel(piezoms1)
        leg{ip}=sprintf('r=%d',piezoms1(ip).vertex(piezoms1(ip).iColx)); %#ok
        piezoms1(ip).plotHead(ax,gr,Hds,[piezoms1(ip).vertexTxt{iclr} piezoms1(ip).vertexTxt{ipat}],'lineWidth',piezoms1(ip).vertex(iwid));
    end
    legend(leg{1:4},3); % only use the first 5 legends (distance info)
end

%% Compute the sum of all in and output terms except wells
labels    = B(1).label; for j=1:numel(labels), labels{j}(labels{j}==' ')=''; end
newLabels = {'STORAGE','DRAINS','ET'};
%newLabels = {'STORAGE','DRAINS','ET','HEADDEP','RECHARGE'};
I = strmatchi(newLabels,labels);
QdiffA = Qdiff(:,I,:);
squeeze(sum(QdiffA,2))   %show that sum over newLabels equals the extraction of all three models

altLabels={'BERGING','VERMINDERDE AFVOER','VERMIDERING VERDAMPING'};

yLim = [-100000 100000];

figure('name',sprintf('water balance, scenario %d',scenario),'pos',screenPos(0.5));
Npl = numel(newLabels)+1;
for ia=Npl:-1:1
    ax(ia) = subplot(Npl,1,ia,'nextPlot','add','xGrid','on','yGrid','on');
    if ia==1
        area(ax(ia), timeYr,sum(QdiffA(:,:,1),2),'faceColor',mf_color(ia));
        title(ax(ia),sprintf('ONTTREKKING %.1f Mm3/a, (scenario %d)',aYear * abs(mean(well(1).Q))/1e6,scenario));
    else
        area(ax(ia), timeYr,QdiffA(:,ia-1,1),'faceColor',mf_color(ia));
        title(ax(ia),altLabels{ia-1});
    end
    set(ax(ia),'yLim',yLim);
    ylabel(ax(ia),'m3/d'); xlabel('tijd [jr]');
end

print(gcf,sprintf('terwisscha_balance_scen_%d.png',scenario),'-dpng'); %#ok


QdiffB = NaN(size(QdiffA(1:10,:,:)));
for i=10:-1:1
    QdiffB(i,:,:) = mean(QdiffA((i-1)*26+1:i*26,:,:),1);
end

for k=1:3
    fprintf([repmat(' %10.0f',[1 size(QdiffB,2)]),'\n'],squeeze(Qdiff(:,:,k))');
end
%%

% show that the sum for each of the three models matches the extraction
Balance = squeeze(sum(QdiffA,2));  % must equal Q

if 0
    zoneArray = zeros(gr.size);
    for i=1:3, zoneArray((2*i-1)+[0 1],:,:)=i; end 
    %for i=1:6, zoneArray(i,:,:)=i; end 
    %zonebudget(B,zoneArray,[1]);
else
    zoneArray = zeros(size(DB(1).term{1}));
    for i=1:3, zoneArray(i,:,:)=i; end 
    %for i=1:6, zoneArray(i,:,:)=i; end 
    %zonebudget(DB,zoneArray,[1]);
end

%% Analyzing linearity of drawdown by comparing dd of Q=2.5 and 7.4Mm3/a

if movie3

    leg = {'Layer 1','Layer 2'};

    if gr.Nlay>2, leg = [leg {'Layer3'}]; end

    time = [Hds.time];  

    figure('pos',screenPos(0.75));
    ax1 = axes('nextplot','add','ylim',[-0.2 1],...
        'xlim',[1000 R1]/1000,'yscale','lin','xscale','log','yDir','reverse',...
        'yGrid','on','xGrid','on','yMinorGrid','off','fontSize',fsz,'xtick',1:R1km);
    ax2 = axes('nextplot','add','ylim',[2 5],...
        'xlim',[1000 R1]/1000,'yscale','lin','xscale','log','yDir','reverse',...
        'yGrid','off','xGrid','off','yMinorGrid','off','fontSize',fsz,'xtick',1:R1km,...
        'color','none','yAxisLocation','right');

    axWeather = axes('pos',sunPos,'nextPlot','add','color','none','visible','off','yDir','reverse');
    Sun  = imread('Sun.png');
    Rain = imread('Rain','png');
    hSun = image(Sun,'parent',axWeather);
    hRain= image(Rain,'parent',axWeather);

    Q1 = abs(mean(well(1).Q))*aYear/1e6;
    Q2 = abs(mean(well(3).Q))*aYear/1e6;
    
    vidObj = VideoWriter([basename '25+75' sprintf('_scen_%d',scenario)]);
    vidObj.FrameRate = 15;
    vidObj.Quality   = 80;
    vidObj.open();
    for it=1:numel(time); % (1:end-26))
        s = sprintf('Drawdown Terwisscha, scenario %d, XSection, Q=%.1fMm3/a vs Q=%.1fMm3/a t=%.2f y',scenario,Q1,Q2,timeYr(it));
        if it==1
             ht1  = title(ax1,s);

             hp7 = plot(ax1,gr.xkm,Dd(it).values(3,:,2),'b');
             hp5 = plot(ax1,gr.xkm,Dd(it).values(3,:,1),'r');
             if gr.Nlay>2
                 hp9 = plot(ax1,gr.xkm,Dd(it).values(3,:,3),'k');
             end
             hp6 = plot(ax1,gr.xkm,Dd(it).values(1,:,1),'r');
             hp8 = plot(ax1,gr.xkm,Dd(it).values(1,:,2),'b');
             if gr.Nlay>2
                 hp10= plot(ax1,gr.xkm,Dd(it).values(1,:,3),'k');
             end
             if gr.Nlay>2
                 hp13 = plot(ax2,gr.xkm,Dd(it).values(1,:,3)./Dd(it).values(3,:,3),'k--','lineWidth',2);
             end
             hp12 = plot(ax2,gr.xkm,Dd(it).values(1,:,2)./Dd(it).values(3,:,2),'b--','lineWidth',2);
             hp11 = plot(ax2,gr.xkm,Dd(it).values(1,:,1)./Dd(it).values(3,:,1),'r--','lineWidth',2);

             plot(ax1,[R,gr.xGr(end)]/1000,[0 0],'g','linewidth',2);
             text(Rkm,0,'Agricultural land -->','VerticalAlignment','Bottom','parent',ax1,'fontsize',fsz);

             legend(ax1,leg,2);

 %            gr.plotXSec(ax1,1,'all','hlines','k','edgecolor',grey,'faceColor','none');   

             xlabel(ax1,'r [km]');
             ylabel(ax1,'drawdown [m]');
             ylabel(ax2,'drawdown7.5Mm3pj/drawdown2.5Mm3pj')
             text(1.1,3,'Ratio','fontSize',16,'parent',ax2,...
                 'fontWeight','bold','fontAngle','italic','verticalAlignment','Bottom');             

        else
            set(ht1,'string',s);
            set(hp5,'ydata',Dd(it).values(3,:,1));  % Q=7.5Mm3/a
            set(hp6,'ydata',Dd(it).values(1,:,1));  % Q=2.5Mm3/a
            set(hp7,'ydata',Dd(it).values(3,:,2));  
            set(hp8,'ydata',Dd(it).values(1,:,2));
            if gr.Nlay>2
                set(hp9,'ydata',Dd(it).values(3,:,3));  
                set(hp10,'ydata',Dd(it).values(1,:,3)); 
            end
            set(hp11,'ydata',Dd(it).values(1,:,1)./Dd(it).values(3,:,1));
            set(hp12,'ydata',Dd(it).values(1,:,2)./Dd(it).values(3,:,2));
            if gr.Nlay>2
                set(hp13,'ydata',Dd(it).values(1,:,3)./Dd(it).values(3,:,3));
            end
        end
        if sunny(it);
            set(hSun,'visible','on'); set(hRain,'visible','off');
        else
            set(hSun,'visible','off'); set(hRain,'visible','on');
        end
        
        vidObj.writeVideo(getframe(gcf));
    end

    vidObj.close();
end

%% Time series analysis

tsanalysis;

%% Menyanthes
% [~,tpe] = getPeriods(basename,'PER',{'IPER','RECH','EVTR'},'fill',true);
% H = piezoms1.menyanthesH(gr,Hds);
% IN = menyanthesIn('tpe',tpe);
% 
% xP = 201402;
% yP = 564795;
% 
% for i=numel(H):-1:1
%     H(i).xcoord = -H(i).xcoord + xP;
%     H(i).ycoord =  H(i).ycoord + yP;
%     H(i).values(:,1) = H(i).values(:,1) + datenum(1987,4,21);
% end
% for i=numel(IN):-1:1
%     IN(i).xcoord = -IN(i).xcoord + xP + 1000;
%     IN(i).ycoord =  IN(i).ycoord + yP + 1000;
%     IN(i).values(:,1) = IN(i).values(:,1)+datenum(1987,4,21);
% end
% 
%     VER = { { '2'    'x'    'e.b'    '(beta)'} '5 December 2013' };
% 
%     T = struct('struct',{},'location',{},'data',{},'date',{},'deleted_by',{});
% 
%     B = struct(...
%         'ID',{},'name',{},'xcoord',[],'ycoord',[],'surflev',{},'date',{},'soillayers',{});
% 
% 
%     ID = struct('H',{},'IN',{},'M',{},'B',{});
% 
%    
%      M = struct(...
%         'ID',{},'name',{},'tdiscr',{},'equid',{},'nonlin',{},'resamph',{},'startup',{},...
%         'lockpar_yn',{},'lockedpar',{},'xcoord',{},'ycoord',{},'zcoord',{},...
%         'h',[],'result',{},'noise',{},'in',{});
% 
% %     IN = struct(...
% %         'ID',{},'name',{},'type',{},'values',{},'xcoord',{},'ycoord',{},...
% %         'surflev',{},'filtnr',{},'upfiltlev',{},'lowfiltlev',{},'datlog_serial',{});
% 
%     
%        ID = struct('H',numel(H),'IN',numel(IN),'M',numel(M),'B',numel(B));
% 
%        IN= [IN IN];
%        IN(end-1).name = 'PREC01_zero';
%        IN(end  ).name = 'EVAP01_zero';
%        IN(end-1).values(:,2)=mean(IN(1).values(:,2));
%        IN(end  ).values(:,2)=mean(IN(2).values(:,2));
%        
% save MyData B H ID IN M T VER
% 
