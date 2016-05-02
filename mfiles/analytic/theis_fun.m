function [s, sac]=theis_fun(Q,K,H,S,r,t    )
% function [s, sac]=theis_fun(Q,K,H,S,r,t    )
% we decide to preserve K H rather than T for the purpose to show all the 
%   processes
% s is approximated solution
% sac is accurate solution;
% theis function
% for sac 
% sac( time index, space index)
% Q --   pumping rate (m3/s
%
%Q_m3_hr=80;
%Q =Q_m3_hr/3600;      %m3/s   discharge rate 80 m3/hr * 1/3600 hr/s
%K =   2.5e-5;  % hydraulic conductivity (m/s)
%H =   200 ;    % thickness of aquifer
%T =K*H;  % transmmisivity     0.00001m/s = 0.864 m/day
%
%S =0.0001;  % storativity
%r = [0.1:0.2:1]; % (m) distance from puming well to the point where 
%                 % draw down was observed
%		 t = [0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 ]*3600*24 ;
%		                  % (seconds) time




s   = struct;    % drawdown
T   = K*H;
s.Q = Q;
s.T = T;
s.K = K;
s.H = H;
s.S = S;
s.t = t;
s.r = r;
s.ttl=['K' num2str(s.K*3600*24) 'm_D__Q' num2str(s.Q*3600) ...
      'm3_h__H' num2str(s.H) 'm__S' num2str(s.S) ]  ;

for i =1:length(t)   % time loop
    u         =  r.^2 * S /4 /T  /t(i);  % at the moment u is an array

    for j=1:length(u) % distance loop j
       summ   = -0.577216-log( u(j) );
       summP  = inf;
       n      = 1;

       while abs(summP-summ)>1e-16
         summP  = summ;
         summ   = summ+ (-1)^ (n-1) * u(j)^n / n / factorial(n);
         n      = n + 1;
       end   % while loop

       s.nloop(i,j)=n;   % number of loops per simulation
       s.dd(i,j)=summ  *  Q  /4 /pi / T  ;
    end  % distance loop j
end % time loop i

sac      = struct;  % accurate result of theis equation;
sac.Q    = Q;
sac.T    = T;
sac.K    = K;
sac.H    = H;
sac.S    = S;
sac.t    = t;
sac.r    = r ;
sac.ttl=['K' num2str(sac.K*3600*24) 'm_D__Q' num2str(sac.Q*3600) ...
      'm3_h__H' num2str(sac.H) 'm__S' num2str(sac.S) ]  ;
for i  = 1: length(t)
    sac.dd(i,:) = Q/ (4*pi*T ) * expint( r.^2 * S /  (4*T*t(i) )  );
    % a wrong solution by Ling see 150731 in book ejicamp
    %sac.dd(i,:) = Q/ (4*pi*T ) * ei( r.^2 * S /  (4*T*t(i) )  );
end

