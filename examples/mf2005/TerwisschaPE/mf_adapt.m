% CIE 5440 23 April 2013, Terwisscha asignment
% Adapted to prepare the presentation for NHV 2013-11-27 in Utrecht
% regarding the article by Van den Akker, Stromingen (19)2.

% For this simple case we use a two layer model that is axially symmetric.
% The top layer is the cover layer and the second is the regional aquifer.

% Drains will be used as areal boundaries because they will not allow
% infiltration when the water table sinks to below their elevation.

basename = 'terwisscha';
save name basename

%% Parameters
ayear = 365.24;     % [d/y]

kD1 =   50; D1 =  30; c1 =  100;  % cover layerlayer
kD2 = 5100; D2 = 100; c2 = 2000;  % first aquifer
kD3 = 5100; D3 = 100;            % second aquifer

% scenarios
vdAkker_National   = false;
vdAkker_Terwisscha = false;

% Defaults
        extDepth = 1.5;  % [m] extension depth in PER sheet
        R  = 3000;   % [m] radius of nature area
        R1 = 8500;
        Nlay=2;
        syNature= 0.25;         % [-] specifid yield in nature area
        syAgri  = 0.15;         % [-] specific yield in aricultural area

[PERnams,PERvals] = getPeriods(basename,'PER');
exDp = mean(PERvals(:,strmatchi('EXDP',PERnams)));
evtr = PERvals(:,strmatchi('EVTR',PERnams));
rech = PERvals(:,strmatchi('RECH',PERnams));
sun  = evtr-rech>0;

scenario = 2;
        
% Special, only what differs from defaults
switch scenario
    case 1 % vdAkker, De Glee, Dupuit Theis
        R  = 1;        
    case 2 % Defaults
    case 3 % Sy in nature and agricltural land
        % make sure extension depth is set to 150 m instead of 1.5 m
        extDepth = 150;
        syNature= 0.10;
        syAgri  = 0.10;
    case 4 % Deep aquifer included
        Nlay = 3;
    case 5 % Model radius at 20000
        R1 = 20000;
    case 6 % Model radius at 2000, deep aquifer included
        Nlay=3;
        R1 = 20000;
end

% Make sure that ET reduction is set correctly for case
if extDepth~=exDp
    error('Change EXDP in PER sheet from %g to %g for scenario %d',exDp,extDepth,scenario);
end

hk      = [kD1/D1  kD2/D2  kD3/D3]; % [m/d]  horizontal conductivit of cover layer

ss      = 1e-5;         % [1/m]; specific elastic storage coeffiient of all layers

Omega   =  1.5;        %  [m] wet circumference of ditches
w       =  0.5;         % [d] entry resistance of ditches

b       =   70;         % [m] half distance between parallel ditches

iWell  =   5;           % [-] zone number for well
iNature=   6;           % [-] zone number nature area
iAgri  =   7;           % [-] zone number agricultural area

% Entry resistance of agricultural land, taking into account both ditches
% and their entry resistance, as is documented in the assignment.
lambda1 = sqrt(kD1*c1);

% Resistane between ditches and phreatic aquifer
cAgri   = b*(w/(Omega/2) + log(pi*D1/Omega)/(pi*hk(1))) + c1*((b/lambda1)*cosh(b/lambda1)/sinh(b/lambda1)-1);

%% Grid
Top= 10;

%rGr = [  0 logspace(0,log10(R1),121)];  % grid from 1 to 10000 m
rGr = sinespace(0,R1,121,0,pi/2);
% Top of model is >0, to allow for the free water surface
% the starting heads 0
% So the average thickness of the top layer is still d.
% Make sure LAYCON==1 in sheet LAY to allow for water table behavior
zGr = Top - cumsum([0 D1 D2 1 D3]) ;      % layer elevations

if Nlay<=2
    % in case of two layer only zGr(1:3) is used. That is the top layer is
    % considered sits directly on the top of the regional aquifer.
    % It also implies that the drainage resistance is combined with the
    % vertical resistance of the cover layer.
    gr = gridObj(rGr,[-5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5],zGr(1:Nlay+1),'AXIAL',true); % well (y=0) is in first row
else
    % In this case we add a 1 m thick aquitard and a second regional
    % aquifer. Hence zGr = zGr(1:5)
    gr = gridObj(rGr,[-5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5],zGr(1:Nlay+2),'AXIAL',true,'LAYCBD',[0 1]); % well (y=0) is in first row
end

if vdAkker_National==true
    if Nlay<=2
        gr = gridObj([0 pi*R1^2],gr.yGr,zGr(1:Nlay+1),'LAYCBD',gr.LAYCBD); % not axial
    else
        gr = gridObj([0 pi*R1^2],gr.yGr,zGr(1:Nlay+2),'LAYCBD',[0 1]); % not axial
    end
end

% Boolean to indicate nature and agricultural area
agri       =  gr.Xm>=R;
nature     = ~agri;
drnMdl     = ismember(gr.Ym,[0 -1 -4 -5]);
ghbMdl     = ~drnMdl;
%% Model arrays

IBOUND     = gr.const(ones(Nlay,1)); % no fixed head boundaries
HK         = gr.const(hk(1:Nlay));

% Only the resistance of the top laer is built into VK not as VKCB
VK         = gr.const(hk(1:Nlay)); VK(:,:,1) = D1/2/c1;

% The aquitard between second and third model layer is implemented as VKCB
VKCB       = gr.constCB(1/c2);
SS         = gr.const(ones(Nlay,1)*ss);
SY         = gr.const(zeros(Nlay,1)); SY(agri) = syAgri; SY(nature) = syNature;
STRTHD     = gr.const(zeros(Nlay,1));

well = wellObj(basename,'wells',gr,HK,{'PER','Q'},'fill',true);

%%
% Boundary conditions for drains
% The conductance CDr of a drain in a cell if the drain is uniformly
% distributed equals CDr = gr.Area/c2;
CDr = NaN(size(gr.AREA));
CDr(agri(  :,:,1)) = gr.AREA(agri(  :,:,1))/cAgri;

zoneValsDRN = {true,  0, CDr(agri & drnMdl)};
zoneValsGHB = {true,  0, CDr(agri & ghbMdl)};

DRN = bcnZone(basename,'DRN',agri & drnMdl,zoneValsDRN);
GHB = bcnZone(basename,'GHB',agri & ghbMdl,zoneValsGHB);

%% De Glee

sDeGlee = well(1).Q(1)/(2*pi*(kD1+kD2))*besselk(0,gr.xm/sqrt(kD2*cAgri));

%% Apply Van den Akker(2013) numerically using fdm2c
%  Figure 9 of Olsthoorn (2014) Tussen De Glee en Dupuit revisited.

vdAkker_Terwisscha=false;
if vdAkker_Terwisscha
    c1=100;
    VK         = gr.const(hk(1:Nlay)); VK(:,:,1) = D1/2/c1;
    FQ = gr.const(0); FQ(1,1,end)=-well(1).Q(1);
    IBND = IBOUND; IBND(agri)=-1;

    cDrainageVdAkker = (D1/2)./(c1+ cAgri * agri); %dummy to rememeber original resistances

    iShow = find(gr.xm<=3500,1,'last');

    figure('name','vgl vdAkker','pos',screenPos(1,0.8));
    xLim = [0.1 9];
    fsz = 16;
    ax1 = subplot(2,1,1,'nextPlot','add','xLim',xLim,'xScale','log','xGrid','on','yGrid','on','yDir','reverse','yScale','lin','fontSize',fsz);
    ax2 = subplot(2,1,2,'nextPlot','add','xLim',xLim,'xScale','log','xGrid','on','yGrid','on','yDir','reverse','yScale','log','fontSize',fsz);
    title(ax1,{'Drawdown regular and according to vdAkker (2013), Q=6.5Mm3/a';...
           'kD1=50m2/d, kD2=5000m2/d, c1=100d, rAgri=3000m, d0=0.6m';...
           'cas1: cDr=132d, case2: a=0.45m; case3: a=0.1m'});
    title(ax2,'same, but y-axis logarithmic');
    xlabel(ax1,'r [km]'); ylabel(ax1,'drawdown [m]');
    xlabel(ax2,'r [km]'); ylabel(ax2,'drawdown [m]');
    lSt = {'o','-','-'};
    for iLoop=1:3
        d0=0.6;
        switch iLoop
            case 1, a=c1./cAgri*d0;
            case 2, a=c1./cAgri*d0;
            case 3, a=0.1;
        end
        for i=1:10;       
            if iLoop>1 && i>1
                deltaq = Qy./gr.AREA(1,:,1);
                cDrainageVdAkker = c1*d0/a * exp(deltaq./(a./c1));
            else
                cDrainageVdAkker = c1*d0/a;
            end
            fprintf('iLoop=%d, cDrain(r=%.0f)=%.0f,  ',iLoop,gr.xm(iShow),cDrainageVdAkker(min(numel(cDrainageVdAkker),iShow)));
            Ky= VK;
            Ky(:,:,1) = (D1/2)./(c1+bsxfun(@times,cDrainageVdAkker,agri));

            [Phi,Q,Psi,Qx,Qy] = fdm2(gr.xGr,gr.zGr(:),XS(HK(1,:,:)),XS(Ky(1,:,:)),...
                        XS(IBND(1,:,:)),XS(STRTHD(1,:,:)),XS(FQ(1,:,:)),'radial');
            %plot(gr.xm,Phi(end,:),mf_color(i));
        end    
        plot(ax1,gr.xm/1000,Phi(end,:),[mf_color(iLoop) lSt{iLoop}],'lineWidth',2);
        plot(ax2,gr.xm/1000,Phi(end,:),[mf_color(iLoop) lSt{iLoop}],'lineWidth',2);
        fprintf('\n');
%         if iLoop==3
%             fprintf('%12f  %12f\n',[gr.xm;(D1/2)./Ky(1,:,1)-100])
%         end
    end
    legend(ax1,'base case','vdAkker a=0.45m','vdAkker a=0.1m',2);
    legend(ax2,'base case','vdAkker a=0.45m','vdAkker a=0.1m',2);
    set([ax1 ax2],'xTick',[0.1:0.1:0.9 1:1:9]);
end
%%
save underneath b R R1 sDeGlee scenario kD1 kD2 scenario sun