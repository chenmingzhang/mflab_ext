% Interpretation of time-series analysis. The time series generated by the
% groundwater model is fast, as it head increase follows immediately after
% recharge. Also, it has no noise. Therefore, we can model it with an
% impulse response that is an exponential decay.
global PE head

names = {'DRN, Q65', 'DRN, Q00', 'GHB, Q65', 'GHB, Q00'};

[PERnams,PERvals,NPER] =  getPeriods(basename,'PER',{'year','month','day','RECH','EVTR'});
date = datenum(PERvals(:,1:3));

%% Put heads in arrays to easy manipulation

measHds{4} = NaN(gr.Nz,gr.Nx,NPER);
measHds{3} = NaN(gr.Nz,gr.Nx,NPER);
measHds{2} = NaN(gr.Nz,gr.Nx,NPER);
measHds{1} = NaN(gr.Nz,gr.Nx,NPER);

for i=1:numel(measHds)
    for it=1:NPER
        measHds{i}(:,:,it) = XS(Hds(it).values(i,:,:));  % Q=65, DRN
    end
end

%% load meteo

load TPE.txt

PE = PERvals(:,4)-PERvals(:,5);

% add random error to PE
fPE = 0.01; % factor
fHd = 0.01; % m

PE = PE + fPE * abs(PE) .* randn(size(PE));
PE = PE - mean(PE);

%% Initial parameters
ymean =    0;
tau0  =   30;
mu    = 0.10;

par0 = [ymean; log(mu); log(tau0)];

options = optimset('FinDiffType','central');

%% Solve

deep = true;

if deep, layer=2; else layer = 1; end

for i = numel(measHds):-1:1
    fprintf('tsa on measured heads series{%d}\n',i);
    
    par{i}= NaN(3,gr.Nx);    

    for ix=1:gr.Nx  % for all points

        fprintf('heads at x = %.0fm\n',gr.xm(ix));
        
        head = squeeze(measHds{i}(layer,ix,:));
        head = head + fHd * randn(size(head));
        
        par{i}(:,ix) = lsqnonlin(@tsFun,par0,-Inf,Inf,options);
        fprintf('\n');
    end

    % turn log variables back to non-log values
    par{i}(2:end,:) = exp(par{i}(2:end,:));
    
    % add x row on top
    par{i}          = [gr.xm; par{i}];      
end

%% plot results for structured head
figure('pos',screenPos(1,0.8));
axes('nextPlot','add','ylim',[-2 1],'fontSize',fsz);
jys = 2; % row in par of structural level
lst = {'b.-','bs-','r.-','rs-'};
leg = {};
for i=1:numel(measHds)
    leg{i} = names{i}; %#ok
    plot(gr.xm/1000, par{i}(jys,:),lst{i});
    grid on;
end
legend(leg{:});
xlabel('distance r [km]'); ylabel('head [m]');
title(sprintf('Structural level layer %d',layer));

%%
figure('name','dd of structural level','pos',screenPos(1,0.8));

fsz = 16;

yLim = [-0.1 2.5];

% linear scale
ax1 = axes('nextPlot','add','ylim',yLim,'fontSize',fsz,'xGrid','on','yGrid','on','yDir','reverse','xScale','log');
% log scale
leg = {sprintf('DRN, layer %d',layer),sprintf('GHB, layer %d',layer)};

%ax2 = subplot(2,1,2,'nextPlot','add','yLim',yLim,'fontSize',fsz,'xGrid','on','yGrid','on','yDir','reverse','xScale','log');
plot(ax1,gr.xm/1000, abs(par{1}(jys,:) - par{2}(jys,:)),'-b');
plot(ax1,gr.xm/1000, abs(par{3}(jys,:) - par{4}(jys,:)),'-r');
legend(ax1,leg{:},4);

%plot(ax2,gr.xm/1000, abs(par{1}(jys,:) - par{2}(jys,:)),'-b');
%plot(ax2,gr.xm/1000, abs(par{3}(jys,:) - par{4}(jys,:)),'-r');
%legend(ax1,leg{:},4);

xlabel(ax1,'r [km]'); ylabel(ax1,'drawdown [m]'); title(ax1','drawdown structural level');
%xlabel(ax2,'r [km]'); ylabel(ax2,'drawdown [m]'); title(ax2','drawdown structural level');
set(ax2,'xTick',[0.1:0.1:1 2:10]);

plot(ax1,gr.xm/1000,drawd1Myr(2,:),'k.');
%plot(ax1,gr.xm/1000,drawd2myr(2,:),'kv');
%plot(ax2,gr.xm/1000,drawd1Myr(2,:),'k.');
%plot(ax2,gr.xm/1000,drawd2myr(2,:),'kv');

%leg=[leg,'minDrawdown','maxDrawdown'];
leg=[leg,'maxDrawdown'];
legend(ax1,leg{:},4);
%legend(ax2,leg{:},4);

set(ax1,'xlim',[1 10],'ylim',[-0.1 1],'xTick',1:10);
set(ax1,'yScale','log');
%set(ax2,'yScale','log');

%% dd Maas(2012) fig 30 , digitized from image
ddMaas = [
          85        -2.425
         412        -1.553
         961        -0.947
        1663        -0.793
        1701        -0.843
        1884        -0.772
        1913        -0.656
        2403        -0.581
        2528        -0.648
        3153        -0.390
        3827        -0.287
        4125        -0.241
        4144        -0.349
        4221        -0.366
        5558        -0.374
        5923        -0.291
        7174        -0.199
        7491        -0.303
        8549        -0.366
        8953        -0.328
        9453        -0.336
        9491        -0.224
        ];

    ddAchtergrond = 0.30;

% Maas log/Dupuit curve
ddLogMaas=[
          95        2.375-0.30
        7000        0.004-0.30
    ];

% Maas: s = Q/(2*pi*kD) * log(r/R); --> kD=5200, R=4100;
    
plot(ax1,ddMaas(:,1)/1000,abs(ddMaas(:,2))-ddAchtergrond,'kp','markerFaceColor','k','markerSize',8);
%plot(ax2,ddMaas(:,1)/1000,abs(ddMaas(:,2))-ddAchtergrond,'kp','markerFaceColor','k','markerSize',8);

leg = [leg 'Maas(2012)'];
legend(ax1,leg{:},2);
%legend(ax2,leg{:},2);

%kD = 5200; R=4100;
%plot(ax1,gr.xm(gr.xm<=R)/1000,well(1).Q(1)/(2*pi*kD)*log(gr.xm(gr.xm<=R)/R),'k','linewidth',2);
%plot(ax2,gr.xm(gr.xm<=R)/1000,well(1).Q(1)/(2*pi*kD)*log(gr.xm(gr.xm<=R)/R),'k','linewidth',2);

%leg = [leg 'Maas(2012)-Dupuit'];
%legend(ax1,leg{:},2);
%legend(ax2,leg{:},2);

set([ax1 ax2],'yScale','lin')