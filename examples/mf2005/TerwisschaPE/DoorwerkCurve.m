%% Doorwerkcurve Van den Akker
% grafiek van de doorwerkcurve afhankelijk van initiele grondwaterdiepte en
% parameter a.
% TO140102.

d = 0:0.025:2;
a = 0.10;

F=@(d,a) d./(a+d);


figure; axes('nextPlot','add','yDir','reverse');
a= 0.1; plot(F(d,a),d,'b','lineWidth',3);

a= 0.6; plot(F(d,a),d,'r','lineWidth',3);


%% Van den Akker

% h = -2:0.02:0; mv=0; b=0;
% A = @(h,mv,b,a) 1./(1+a./(h-mv-b)); 
% 
% a=0.1; plot(A(h,mv,b,-a),mv-h,'go');
%a=1.0; plot(A(h,mv,b,-a),mv-h,'gx');

%% Verlaging op afstand numeriek

kD1 = 
kD2 =
c1  =

[Phi,Q,Psi,Qx,Qy] = fdm2c(gr.xGr,gr.zGr(:),XS(HK(1,:,:)),XS(VK(1,:,:)),...
            XS(IBOUND(1,:,:)),XS(STRTHD(1,:,:)),XS(FQ(1,:,:)));
        
        
        
