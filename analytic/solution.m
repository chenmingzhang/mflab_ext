classdef solution
    %SOLUTION class definition for objects of anlytical solutions for dynamics in 1D top layer
    %
    % Class def for anlytical solution objects to compute dynamics of grw fluct. in a 1 and 2 layer aquifer
    % Inlcludes one numeric solution for the dynamics of the average head in a one- or twolayer cross section
    % TO 110714 110917 (now works for single cross section)
    properties
        name      % name of the solution
        b         % [ m ] half width of parcel
        dd        % [ m ] ditch depth
        wd        % [ m ] ditch width
        GR_ELEV   % [NAP] elevation of ground surface
        t         % [ d ]time
        N         % [m/d] net recharge
        q         % [m/d] seepage (upward positive)
        D         % [ m ] effective thickness of the two aquifers
        Kx        % [m/d] horizontal conductivity in the two aquifers
        Ky        % [m/d] vertical   conductivity in the two aquifers
        kD        % [m2/d] transmissivity of the two aquifers
        c         % [ d ]  resistance between the two aquifers
        lambda    % [ m ]  sqrt(kD(:,1).*c)
        gamma     % [ - ]  b/lambda
        Lambda    % [ - ]  tanh(gamma)/gamma
        Gamma     % [ - ]  Gamma/(1-Gamma)
        S         % [ - ]  storage coefficient of the two aquifers
        win       % [ t ]  entry resistance for infiltrating ditch (only implemented for solution semiw (experimental)
        wex       % [ t ]  exfiltration resistance of ditch (check definition), exfiltration=infiltration=default
        cex       % [ t ]  exfiltration resistance of ditch bottom and side
        cin       % [ t ]  infiltration resistance of ditch bottom and side
        h_mean    % [NAP]  average ditch level
        h_summer  % [NAP]  summer  ditch level
        h_winter  % [NAP]  winter  ditch level
        Qsw       % surface water runoff
        ht        % transient mean head
        hx        % steady head on input of last time step
        he        % mean final head based on last time step
        xGr       % [ m ] cell boundary coordinates along x-axis
        xm        % [ m ] cell center coordinates along x-axis
        Z         % [NAP] top and bottom of the 2 aquifes (2 values per section)
        HY        % GXG struct (contains GXG data for each year in time series
        GLG       % overall GLG
        GVG       % overall GVG
        GHG       % overall GHG
    end
    properties (Constant)
        solutions={
            'Conv'   'Confined exact by convolution';
            'Cos'    'Confined cos approx';
            'Conf'   'Confined seepage added to recharge';
            'Semi'   'Semiconf on aquitard on constant head specified average seepage';
            'GGOR'   'Semiconf on aquitard on constant head specified average seepage, GGOR Ouboter';
            'Semiw'  'Semiconf on aquitard on constant head specified seepage entry resistance;'                    
            'Numeric'      '2 Layer numerical fdm with entry resistances';
            'Semi2w_ode45' '2 layer solution using numerical ode45 solver';
            'Semi2w_rk3'   '2 layer solution using numericla 3rd order Runge-Kutta method';
            'Semi2w_xi'    '2 layer solution direct integration using intermediate variable xi';
            'Semi2w'       '2 layer solution direct integration';
            };
    end
    methods
        function obj=solution(P)
            if nargin==0, obj.name='empty'; return; end
            if length(P)>1,
                obj(length(P),1)=solution();
                for i=1:length(P)
                  obj(i)=solution(P(i));
                end
                return;
            else
                % [  ]; [ ]]' constructions were meant to work for multiple
                % cross sections, but are no longer in use. They work also
                % for single cross sections
                obj.name=    '';
                obj.b   =    P.L/2;
                obj.GR_ELEV =P.AHN;
                obj.t   =    [];
                obj.N   =    [];
                obj.q   =    P.q;
                %obj.D   =    [[P.D1]-([P.AHN]-[P.GP]); P.D2]';
                obj.D   =    [P.D1effective  P.D2];
                obj.Kx  =    [P.hk1 P.hk2];
                obj.kD  =    obj.Kx.*obj.D;
                obj.c   =    P.C;
                obj.Ky  =    [0.5*obj.D(1)./obj.c 0.2*obj.Kx(2)];
                obj.dd  =    P.dd;  % ditch depth
                obj.wd  =    P.wd;  % ditch width
                
                obj.lambda = sqrt(obj.kD(:,1).*obj.c);
                obj.gamma  = obj.b./obj.lambda;
                obj.Lambda = tanh(obj.gamma)./(obj.gamma);
                obj.Gamma  = (3./obj.gamma.^2).*(1-obj.Lambda)./obj.Lambda;

                obj.S   =    [P.MU P.S2];
                obj.cex =    P.cex;
                obj.cin =    P.cin;
                
                % see theory for background
                obj.wex =    [obj.cex*obj.D(1)/min(obj.dd+obj.wd/2,obj.D(1))
                              obj.D(2)/pi/sqrt(obj.Ky(2)*obj.Kx(2))*log(obj.D(2)/obj.wd*(obj.Ky(2)/obj.Kx(2))^2)+ ...
                                (obj.cex+max(0.01,obj.D(1)-P.dd)/obj.Ky(1))*obj.D(2)/obj.wd/2];
                % see theory for background
                obj.win =    [obj.cin*obj.D(1)/min(obj.dd+obj.wd/2,obj.D(1))
                              obj.D(2)/pi/sqrt(obj.Ky(2)*obj.Kx(2))*log(obj.D(2)/obj.wd*(obj.Ky(2)/obj.Kx(2))^2)+ ...
                                (obj.cin+max(0.01,obj.D(1)-obj.dd)/obj.Ky(1))*obj.D(2)/obj.wd/1];
                            
                obj.h_mean  = P.GP;
                obj.h_summer= P.ZP;
                obj.h_winter= P.WP;
                obj.Qsw = [];   % surface water runoff
                obj.ht  = [];   % transient mean head
                obj.hx  = [];   % steady head on input of last time step
                obj.he  = [];   % mean final head based on last time step
                obj.xGr = [];   % coordinates of x-axis
                obj.xm  = [];
                obj.Z   = [];
                obj.HY  = [];   % GVG data struct
                obj.GLG = [];   % overall GLG
                obj.GVG = [];   % overall GVG
                obj.GHG = [];   % overall GHG

                z=repmat(obj.h_mean,[1,1,3]);
                z(:,:,2) = z(:,:,1)-obj.D(:,1);
                z(:,:,3) = z(:,:,2)-obj.D(:,2);
                obj.Z=z;
            end
        end
        function obj=solve(obj,name,tne,xGr)
            obj=obj.legal(name);
            
            % Grid here, because we need xGr also for the analytical
            % steady state solution, not only for the numerical one
            if ~exist('xGr','var')
                dx=1;    xGr=[0.001 0:dx:max(obj.b) obj.b];
            end
            [obj.xGr,~,obj.xm,~,Dx,~,Nx,~]=modelsize(xGr,1:3);           
            Nz = size(obj.Z,3)-1;
                        
            Dt=[NaN; diff(tne(:,1))];  Dt(1)=Dt(2);
            obj.N=tne(:,2)-tne(:,3);  % recharge
            Nt=length(Dt);
            obj.Qsw=zeros(size(obj.b),Nt);

            % summer and winter ditch level
            obj.t=tne(:,1);
            [~,MON]=datevec(tne(:,1));
            
            obj.ht=zeros(Nz,length(tne(:,1)));
            obj.hx=zeros(Nz,length(obj.xm));
            obj.he=zeros(Nz,1);
                      
            switch obj.name
                case 'Conv'
                    % Single aquifer confined by convolution derived from
                    % Kraaijnhof van der Leur, series solution. Exact.
                    fprintf('Warning, Conv cannot compute surface runoff!\n');

                    T    =4*obj.b^2*obj.S(1)/(pi^2*obj.kD(1));

                    n=14; tau  =0:Dt(1):n*T; % derived n=7, taken 2x to be sure

                    hLR = ones(size(tne(:,1))) * obj.h_summer; hLR(tne(:,end)==0)= obj.h_winter;
                    
                    % here we do time at once, but loop over the cross
                    % sections
                    obj.ht=hLR(1)+...
                        filter(BR_recharge(T,obj.S(1),tau),1,obj.N+obj.q)'+...
                        filter(BR_ditch(T,tau),1,hLR(:)-hLR(1))';
                    
                    % steady state
                    F=hLR(end)+3/2*(obj.N(end)+obj.q).*(obj.b.^2./(3*obj.kD(:,1)));
                    for ix=1:length(obj.xm)
                        obj.hx(:,ix)=F.*(1-(obj.xm(ix)./obj.b).^2);
                    end
                    obj.he =((obj.N(end)+obj.q).*obj.b^2)./(3*obj.kD(:,1)); 

                    
                case 'Cos'  % Single aquifer cosine approx, i.e. first term of
                            % Kraaijenhof van der Leur. Seepage directly added to recharge
                    
                    T    =4/(pi^2)*obj.b.^2.*obj.S(1)./obj.kD(1);

                    for it=1:length(Dt);
                        hLR=obj.h_winter+(MON(it)>=4 && MON(it)<=9).*(obj.h_summer-obj.h_winter);
                        e=exp(-Dt(it)./T);
                        if it==1
                            obj.ht(:,it)=hLR + (obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        else
                            obj.ht(:,it)=hLR +(obj.ht(:,it-1)-hLR).*e + (obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        end
                        % tackle surface water runoff
                        if obj.ht(1,it)>obj.GR_ELEV
                            obj.Qsw(it)=obj.b.*obj.S.*(obj.ht(1,it)-obj.GR_ELEV);
                            obj.ht(1,it)=obj.GR_ELEV;
                        end
                    end
                    
                case 'Conf'  % Single confined aquifer seepage directly added to recharge.
                    
                    T   = obj.b.^2.*obj.S(1)./(3*obj.kD(1));  % characteristic time of this solution
                    
                    for it=1:length(Dt);
                        hLR=obj.h_winter+(MON(it)>=4 && MON(it)<=9).*(obj.h_summer-obj.h_winter);
                        e=exp(-Dt(it)./T);
                        if it==1
                            obj.ht(:,it)=hLR+(obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        else
                            obj.ht(:,it)=hLR +(obj.ht(:,it-1)-hLR).*e + (obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        end
                        % tackle surface water runoff
                        if obj.ht(1,it)>obj.GR_ELEV
                            obj.Qsw(it)=obj.b.*obj.S(1).*(obj.ht(1,it)-obj.GR_ELEV);
                            obj.ht(1,it)=obj.GR_ELEV;
                        end
                    end
                    
                case 'Semi'  % Aquif on top of aquitard below which constant head

                    T      =obj.b.^2.*obj.S(1)./(3*obj.kD(1)).*obj.Gamma;
                    
                    for it=1:length(Dt);
                        hLR=obj.h_winter+(MON(it)>=4 && MON(it)<=9).*(obj.h_summer-obj.h_winter);
                        e=exp(-Dt(it)./T);
                        if it==1
                            obj.ht(:,it) =hLR +(obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        else
                            obj.ht(:,it)=hLR +(obj.ht(:,it-1)-hLR).*e + (obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        end
                        % tackle surface water runoff
                        if obj.ht(1,it)>obj.GR_ELEV
                            obj.Qsw(it)=obj.b.*obj.S(1).*(obj.ht(1,it)-obj.GR_ELEV);
                            obj.ht(1,it)=obj.GR_ELEV;
                        end
                    end
                    
                    % steady state
                    R=(obj.N(end)+obj.q).*obj.c;
                    P=cosh(obj.gamma);
                    for ix=1:length(obj.xm)
                        obj.hx(ix)   =R.*(1-cosh(obj.xm(ix)./obj.lambda)./P);
                    end
                    obj.he   =hLR+R.*obj.Gamma;
                    
                case 'GGOR' % single aquifer on top of aquitard below which constant head
                            % As applied in the GGOR-tool bij Waternet (till 2011)
                            % TO 101113 101114 110106

                    EPS=0.9; EPS_m1=1/EPS; % implicitness
                    F = obj.Lambda./(1-obj.Lambda);
                    
                    for it=1:length(Dt)
                        hLR=obj.h_winter+(MON(it)>=4 && MON(it)<=9).*(obj.h_summer-obj.h_winter);
                        Fmuc=obj.c.*obj.S(:,1)/Dt(it);
                        if it==1,
                            obj.ht(:,it)=EPS_m1*(hLR.*Fmuc+hLR.*F+...
                                obj.c.*(obj.N(it)+obj.q))./(F+Fmuc)+(1-EPS_m1)*hLR;
                        else
                            obj.ht(:,it)=EPS_m1*(obj.ht(:,it-1).*Fmuc+hLR.*F+...
                                obj.c.*(obj.N(it)+obj.q) )./(F+Fmuc)+...
                                    (1-EPS_m1)*obj.ht(:,it-1);
                        end
                        % tackle surface water runoff
                        if obj.ht(1,it)>obj.GR_ELEV
                            obj.Qsw(it)=obj.b.*obj.S(1).*(obj.ht(1,it)-obj.GR_ELEV);
                            obj.ht(1,it)=obj.GR_ELEV;
                        end                       
                    end
                    
                    % steady state
                    R=(obj.N(end)+obj.q).*obj.c;
                    P=cosh(obj.gamma);
                    for ix=1:length(obj.xm)
                        obj.hx(ix)   =R.*(1-cosh(obj.xm(ix)./obj.lambda)./P);
                    end
                    obj.he   =hLR+R.*obj.Gamma;

           
                case 'Semiw';                    
                    % Single aquifer on top of aquitad below which constant head
                    % seepage from below specified
                    % The ditches have entry resistance (same as the multi-layer solution)
                    % however, in case the ditch infiltrates a different
                    % resistance is used, yielding B1 and T1 instead of B
                    % and T

                    B = (1 - ...
                        (obj.b./obj.lambda*(obj.Kx(1).*obj.wex(1)./obj.lambda+...
                        1./tanh(obj.b./obj.lambda)))...
                        .^(-1)).^(-1);
                    T = obj.c .* obj.S(1) * 1 ./ (B -1);
                    
                    % in case of infiltration from ditch use different
                    % ditch bottom resistance win
                    B1= (1 - ...
                        (obj.b./obj.lambda*(obj.Kx(1).*obj.win(1)./obj.lambda+...
                        1./tanh(obj.b./obj.lambda)))...
                        .^(-1)).^(-1);
                    T1= obj.c .* obj.S(1) * 1 ./ (B1-1);

                    for it=1:length(Dt) 
                        hLR=obj.h_winter+(MON(it)>=4 && MON(it)<=9).*(obj.h_summer-obj.h_winter);
                        
                        % assume flow from aquifer to ditch
                        e=exp(-Dt(it)./T);
                        if it==1,
                            obj.ht(:, 1)=hLR+(obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e);
                        else
                            obj.ht(:,it)=hLR+(obj.ht(:,it-1)-hLR).*e + (obj.N(it)+obj.q).*T./obj.S(:,1).*(1-e); 
                        end
                        
                        if obj.ht(1,it)<hLR % in case flow from ditch into aquifer
                            e=exp(-Dt(it)./T1);
                            if it==1,
                                obj.ht(:, 1)=hLR+(obj.N(it)+obj.q).*T1./obj.S(:,1).*(1-e);
                            else
                                obj.ht(:,it)=hLR+(obj.ht(:,it-1)-hLR).*e + (obj.N(it)+obj.q).*T1./obj.S(:,1).*(1-e); 
                            end
                        end

                        
                        % tackle surface water runoff
                        if obj.ht(1,it)>obj.GR_ELEV;
                            obj.Qsw(it)=obj.b*obj.S(1).*(obj.ht(1,it)-obj.GR_ELEV);
                            obj.ht(1,it)=obj.GR_ELEV;
                        end
                    end
                    
                   % steady state
                    R=(obj.N(end)+obj.q).*obj.c;
                    P=obj.Kx(:,1).*obj.win(:,1)./obj.lambda.*sinh(obj.gamma)+cosh(obj.gamma);
                    for ix=1:length(obj.xm)
                        obj.hx(:,ix)   =R.*(1-cosh(obj.xm(ix)./obj.lambda)./P);
                    end
                    obj.he   =hLR+R.*T./(B-1);
                    
                case 'Numeric' % Numerical 2Layer solution with entry resistances

                    % problem tackled on a per cross section basis to allow
                    % using the two-D matlab finite difference models

                    obj.ht=zeros(Nz,Nt);                                        
                    obj.hx=zeros(Nz,length(obj.xm));
                    obj.he=zeros(Nz,1);
                    
                    FQ = zeros(Nz,Nx);
                    FH = NaN(Nz,Nx);
 
                    zGr=squeeze(obj.Z(1,:));

                    kx = obj.Kx(:)*ones(size(Dx));
                    C  = obj.c    *ones(size(Dx));
                    ky = Inf;
                    
                    wequiv=obj.win(:)+0.5*Dx(1)./kx(:,1);

                    kx(:,1)=0.5*Dx(1)./wequiv(:);

                    Ss  = (obj.S(:)./obj.D(:))*ones(size(Dx));

                    %% Initial and boundary conditions

                    for it=1:length(Dt)
                        hLR=obj.h_winter+(MON(it)>=4 & MON(it)<=9)*(obj.h_summer-obj.h_winter);
                        FQ(1,:)=obj.N(it).*Dx;
                        FQ(2,:)=obj.q    .*Dx;

                        if it==1
                            IH = ones(Nz,Nx) *hLR;
                            FH(:,1)= IH(:,1);
                            FI=fdm2ct(obj.xGr,zGr,obj.t(1:2),kx,C,ky,Ss,IH,FH,FQ); % compute end of first day separately
                        else
                            FH(:,1)=hLR;
                            FI=fdm2ct(obj.xGr,zGr,obj.t([it-1,it]),kx,C,ky,Ss,FI(:,:,end),FH,FQ);
                        end

                        % tackle surface runoff
                        fi=FI(1,:,end); I=find(fi>obj.GR_ELEV);
                        if ~isempty(I)
                            obj.Qsw(it)=sum((fi(I)-obj.GR_ELEV).*((Ss(1,I)*Obj.D(1)).*Dx(I)));
                            fi(I)=obj.GR_ELEV;                            
                            FI(1,:,end)=fi;
                        end                           
                        obj.ht(1,it)=sum(FI(1,:,end).*Dx,2)/obj.b;
                        obj.ht(2,it)=sum(FI(2,:,end).*Dx,2)/obj.b;
                        
                        % Steady state
                        hx_=squeeze(FI(:,:,end));
                        obj.he(1)=sum(hx_(1,:).*Dx)/obj.b;
                        obj.he(2)=sum(hx_(2,:).*Dx)/obj.b;

                    end
                    fprintf('\n');
                    
                    % simulation results for cross section
                    obj.hx=  fdm2c(obj.xGr,zGr,kx,C,ky,FH,FQ);
                    obj.he = sum(obj.hx.*[Dx; Dx],2)./obj.b;

                otherwise
                    % Analytic solution for single aquifer on top of seminconfined aquifer separated
                    % by leaking semi-confined bed with given supply in second aquifer.
                    %
                    % TO 101113 101114 101211
                    
                    NLay=2; I = eye(NLay);
                    obj.ht=zeros(NLay,length(Dt));
                    obj.hx=zeros(NLay,length(obj.xm));
                    obj.he=zeros(NLay,1);

                        
                        hLR = obj.h_winter+(MON>= 4 & MON<=9)*(obj.h_summer-obj.h_winter);
                        kd=obj.kD(:);
                        cc   = [1e6 obj.c Inf]';
                        s    = obj.S(:);
                        
                        T   =diag(kd); T_m1=T^(-1);
                        s   =diag(s); S_m1=s^(-1);
                        H_m1=diag(obj.Kx(:)'.*obj.win(:)');

                        A=-diag( 1./(kd(2:NLay  ).*cc(2:NLay)),-1)+...
                            +diag( 1./(kd(1:NLay  ).*cc(1:NLay))+1./(kd(1:NLay).*cc(2:NLay+1)), 0)+...
                            -diag( 1./(kd(1:NLay-1).*cc(2:NLay)),1);

                        A_m1  = A^(-1);
                        sqrtA = sqrtm(A); sqrtA_m1=sqrtA^(-1); % sqrtA_m1= sqrtA^(-1);

                        sinhm=(expm(obj.b*sqrtA)-expm(-obj.b*sqrtA))/2;
                        coshm=(expm(obj.b*sqrtA)+expm(-obj.b*sqrtA))/2;

                        F    = (H_m1*sqrtA*sinhm+coshm);
                        F_m1 =F^(-1);
                        TAB    = T*A*(I-sqrtA_m1/obj.b * sinhm * F_m1)^(-1);
                        G=S_m1*TAB;
                        [V,E]=eig(G); V_m1=V^(-1); E_m1=E^(-1);
                        
                        h  =zeros(2,Nt);
                        
                        switch obj.name
                            case 'Semi2w_ode45'
                                fprintf('Warning, Semi2w_ode5 cannot compute surface runoff!\n');

                                [~,h]=ode45('solveNLay',tne(:,1),[hLR(1); hLR(1)],odeset,S_m1,TAB,tne,obj.q,hLR);
                                obj.ht(:)=h';
                                
                            case 'Semi2w_rk3'  %% Diect intration by third order Runge Kutta method with step control  
                                
                                for it=1:length(Dt)
                                    fprintf('.');
                                    hLR = obj.h_winter+(MON(it)>=4 && MON<=9)*(obj.h_summer-obj.h_winter);
                                    if it==1
                                        H0 =[hLR;hLR];
                                    else
                                        H0 = h(:,it-1);
                                    end
                                    dt=Dt(it); tau=0;
                                    while tau<Dt(it)
                                        f1=S_m1*([obj.N(it);obj.q]-TAB*(H0-hLR)); H1=H0+f1*dt/2;
                                        f2=S_m1*([obj.N(it);obj.q]-TAB*(H1-hLR)); H2=H0+dt*(2*f2-f1);
                                        f3=S_m1*([obj.N(it);obj.q]-TAB*(H2-hLR));
                                        H3=H0+dt*(f1+4*f2+f3)/6;
                                        d=dt*(f1-2*f2+f3)/6; 
                                        if  max(abs(d))<1e-4
                                            tau=tau+dt;
                                            H0=H3;
                                        else
                                            dt=dt/2;
                                        end
                                    end
                                    h(:,it)=H3;
                                    % tackle surface water runoff
                                    if h(1,it)>obj.GR_ELEV
                                        obj.Qsw(it)=obj.b*obj.S(1)*(h(1,it)-obj.GR_ELEV);
                                        h(1,it)=obj.GR_ELEV;
                                    end

                                    if rem(it,100)==0, fprintf('%d\n',it); end
                                end
                                
                            case 'Semi2w_xi' % Making the equations independent by means of eigen vectors and eigen values
 
                                xi     =zeros(2,Nt);
                                for it=1:length(Dt)
                                    hLR = obj.h_winter+(MON(it)>=4 && MON(it)<=9)*(obj.h_summer-obj.h_winter);
                                    e=expm(-E*Dt(it));
                                    theta=V_m1*S_m1*[obj.N(it);obj.q];
                                    if it==1
                                        xi(:,it)=(I-e)*E_m1*theta;
                                    else
                                        xi(:,it)=e*xi(:,it-1)+(I-e)*E_m1*theta;
                                    end
                                    h(:,it)=hLR+V*xi(:,it);
                                    % tackle surface water runoff
                                    if h(1,it)>obj.GR_ELEV
                                        obj.Qsw(it)=obj.b*obj.S(1)*(h(1,it)-obj.GR_ELEV);
                                        h(1,it)=obj.GR_ELEV;
                                        xi(:,it)=V_m1*(h(:,it)-hLR);
                                    end
                                end

                            case  'Semi2w'; % omit the step via xi
                                R=E_m1*V_m1*S_m1;
                                
                                for it=1:length(Dt)
                                    hLR = obj.h_winter+(MON(it)>=4 && MON(it)<=9)*(obj.h_summer-obj.h_winter);
                                    e=expm(-E*Dt(it));
                                    if it==1
                                       h(:,it)=hLR+V*(I-e)*R*[obj.N(it);obj.q];
                                    else
                                       h(:,it)=hLR+V*e*V_m1*(h(:,it-1)-hLR)+V*(I-e)*R*[obj.N(it);obj.q];
                                    end
                                    % tackle surface water runoff
                                    if h(1,it)>obj.GR_ELEV
                                        obj.Qsw(it)=obj.b*obj.S(1)*(h(1,it)-obj.GR_ELEV);
                                        h(1,it)=obj.GR_ELEV;
                                    end
                                end
                        end
                        obj.ht(1,:)=h(1,:);
                        obj.ht(2,:)=h(2,:);
                        
                        % steady state
                        hx_=zeros(NLay,length(obj.xm));
                        for ix=1:length(obj.xm);
                            hx_(:,ix)=hLR+(I-funm(obj.xm(ix)*sqrtA,@cosh)*F_m1)*A_m1*T_m1*[obj.N(it);obj.q];
                        end
                        obj.hx=permute(hx_,[3,2,1]);
                        obj.he=sum(hx_.*(ones(Nz,1)*Dx),2)./obj.b;
             end
        end
        
        function obj=getGXG(obj)
            %[GLG,GVG,GHG]=getGXG(h,t,plotmode)
            % Compute Average lowest head (GLG), spring head (GVG) and highest head
            % (GHG) from the given heads over the given time span.
            % the heads are presumed ordere in the direction of the time vector
            % therefore, a large number of time series can be handled at once
            % plot=0 or omitted altogether, don't plot
            % plot=1 plot heads over time with points that were used to compute CxG
            % plot=2 plot GxG for all sections with section number on x-axis
            % TO 101023

            DV=datevec(obj.t(:));
            
            % get hydrological years in time series
            if obj.t(  1)<=datenum(DV(  1,1),4, 1),  % if first date < April 1 in year 1
                yr1=DV(  1,1);
            else
                yr1=DV(  1,1)+1;
            end
            if obj.t(end)>=datenum(DV(end,1),3,31), % if last data > 31 March last year
                yr2=DV(end,1);
            else
                yr2=DV(end,1)-1;
            end

            yrs=yr1:yr2;  % hydrological years in time series

            obj.HY(length(yrs)).t1=NaN; % preallocate memory for hydrological years
            
            for ihy=1:length(obj.HY)    % for all hydrological years in data

                % GVG
                K=find(DV(:,1)==yrs(ihy) & DV(:,2)==4 & DV(:,3)==1); % april first in hydrological years
                obj.HY(ihy).GVG=[obj.t(K),obj.ht(1,K)];

                % GHG and GLG
                obj.HY(ihy).t1=datenum(yrs(ihy)  ,4,1);  % its starting datenum
                obj.HY(ihy).t2=datenum(yrs(ihy)+1,4,1);  % its ending   datenum
                obj.HY(ihy).J=find(obj.t>=obj.HY(ihy).t1 & obj.t<obj.HY(ihy).t2 & (DV(:,3)==14 | DV(:,3)==28));
 
                R=sortrows([obj.t(obj.HY(ihy).J) obj.ht(1,obj.HY(ihy).J)'],2);

                obj.HY(ihy).GLG=R(1:3      ,:);
                obj.HY(ihy).GHG=R(end-2:end,:);
            end
            obj.GHG=mean(vertcat(obj.HY.GHG)); obj.GHG=obj.GHG(2);
            obj.GLG=mean(vertcat(obj.HY.GLG)); obj.GLG=obj.GLG(2);
            obj.GVG=mean(vertcat(obj.HY.GVG)); obj.GVG=obj.GVG(2);
        end
        
        function obj=showGXG(obj)
        %% Show data and points picked out for GXG computation
            yrs=datevec([obj.HY.t1]);

            figure; hold on;
            xlabel(sprintf('%d - %d',yrs(1),yrs(end)));
            ylabel('head [m]'); grid on;
            title(sprintf('Solution=%s: head and values picked out for GXG',obj.name));
            plot(obj.t([1 end]),[obj.GHG obj.GHG],'r--');
            plot(obj.t([1 end]),[obj.GVG obj.GVG],'g--');
            plot(obj.t([1 end]),[obj.GLG obj.GLG],'b--');
            
            J=vertcat(obj.HY.J);
            plot(obj.t(J),obj.ht(1,J), 'bo'); % plot all
            ghg=vertcat(obj.HY.GHG);
            gvg=vertcat(obj.HY.GVG);
            glg=vertcat(obj.HY.GLG);
            
            plot(ghg(:,1),ghg(:,2),'ro','markerFaceColor','r');
            plot(gvg(:,1),gvg(:,2),'go','markerFaceColor','g');
            plot(glg(:,1),glg(:,2),'bo','markerFaceColor','b');
            plot(obj.t,obj.ht(1,:),'k');
            datetick('x',12);

            legend('GHG','GVG','GLG','sim14','GHGyr','GVGyr','GLGyr','sim');
        end

        function obj=plott(obj,clr,linestyle,t)
            twolayers = size(obj.ht,2)>2;
            hold on; grid on;
            if nargin<4, t=obj.t; end
            if nargin<3, linestyle='-'; end
            if nargin<2, clr='r'; end
            plot(t,obj.ht(1,:),[clr(1) linestyle(1)]);
            if twolayers
                plot(t,obj.ht(2,:),[clr(end) linestyle(end)])
            end
            title(sprintf('Solution %s',obj.name));
            xlabel('time'); ylabel('head [m]');
            %datetick;
        end
        
        function obj=plotx(obj,clr)
            twolayers = length(size(obj.hx(:,1)))>1;
            hold on; grid on;
            if nargin==1, clr='br'; end
            plot(obj.xm,obj.hx(1,:),clr(1));
            plot(obj.xGr([1 end]),obj.he(1,[1 1]),[clr(1) '--']);
            if twolayers
                plot(obj.xm,obj.hx(2,:),clr(end));
                plot(obj.xGr([1 end]),obj.he(2,[1 1]),[clr(end) '--']);
            end
            title(sprintf('Solution %s, steady XSec',obj.name));
            xlabel('x [m]'); ylabel('head [m]');
        end
        
        function br=BR_recharge(T,mu,tau)
            % BR_recharge block response for recharge
            NRech=1;
            n=max(30,1.5*T/diff(tau(1:2))); % afgeleide criterium bleek te zwak, 30 door experiementeren
            sr=(1-exp(-tau/T));
            for j=2:n
                j2=2*j-1;
                sr=sr+(1/j2)^4*(1-exp(-j2^2*tau/T));
            end
            sr=NRech*T/mu*8/(pi^2)*sr;
            br=diff(sr);
        end
        
        function br=BR_ditch(T,tau)
            % BR_ditch block response for ditch
            A=1;
            n=max(30,1.5*T/diff(tau(1:2))); % afgeleide criterium bleek te zwak, 30 door experimenteren
            sr=exp(-tau/T);
            for j=2:n
                j2=2*j-1;
                sr=sr+(1/j2)^2*exp(-j2^2*tau/T);
            end
            sr=A*(1-8/(pi^2)*sr);
            br=diff(sr);            
        end
        function obj=legal(obj,name)
            if nargin==1
                fprintf('Legal names for the implemented solutions:\n');
                for j=1:size(obj.solutions,1)
                    fprintf('%15s : %s\n',obj.solutions{j,:});
                end
            else
                I=strcmpi(name,obj.solutions(:,1));
                if any(I)
                    obj.name=obj.solutions{find(I>0,1,'first')};
                else
                    error('Solution name <<%s>> not known\n',name);
                end
            end
        end
    end
end

