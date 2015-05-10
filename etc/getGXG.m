function [GLG,GVG,GHG]=getGXG(h,t,plotmode)
%GETGXG computes average lowest, spring and highest groundwater level
%
% USAGE:
%    [GLG,GVG,GHG]=getGXG(h,t,plotmode)
%
% h=h(NSection,Nt)
%
% Computes average lowest (GLG), spring (GVG) and highest (GHG) ground water level
% from piezometer time series over the given time span.
%
% used in ..
%   mflab/examples/mf2005/DutchTop/HSV
%   mflab/examples/mf2005/DutchTop/NPARK
%
% the heads are presumed ordered in the direction of the time vector.
% Therefore, a large number of time series can be handled at once
% plot=0 or omitted altogether, don't plot
% plot=1 plot heads over time with points that were used to compute CxG
% plot=2 plot GxG for all sections with section number on x-axis
%
% TO 101023

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

DV=datevec(t(:));
if t(  1)<=datenum(DV(  1,1),4, 1), yr1=DV(  1,1); else yr1=DV(  1,1)+1; end
if t(end)>=datenum(DV(end,1),3,31), yr2=DV(end,1); else yr2=DV(end,1)-1; end

yrs=yr1:yr2;

%% Hydrological years 1-4 tot 31-3

HY(length(yrs)).t1=NaN; % preallocate memory
for i=1:length(HY)
    HY(i).t1=datenum(yrs(i)  ,4,1);
    HY(i).t2=datenum(yrs(i)+1,4,1);
    HY(i).I=find(t>=HY(i).t1 & t<HY(i).t2);
    HY(i).J=find(t>=HY(i).t1 & t<HY(i).t2 & (DV(:,3)==14 | DV(:,3)==28));
    HY(i).K=find(DV(:,1)==yrs(i) & DV(:,2)==4 & DV(:,3)==1); % april first in hydrological years

    [h_yr Is]=sort(h(:,HY(i).J),2);   
    HY(i).GLG  = h_yr(:,1:3); % mean of lowest three
    HY(i).GLGJ = Is(:,1:3);
    
    HY(i).GHG  = h_yr(:,end-2:end); % mean of highest three
    HY(i).GHGJ = Is(:,end-2:end);
    
    HY(i).GVG=   h(:,HY(i).K);
    HY(i).GVGJ=  HY(i).K; 
end
GHG=mean([HY.GHG],2);
GLG=mean([HY.GLG],2);
GVG=mean([HY.GVG],2);

%% Show data and points picked out for GXG computation

if size(h,1)==1 && plotmode>0
    figure; hold on; xlabel(sprintf('%d - %d',yrs(1),yrs(end))); ylabel('head [m]'); grid on;
    title('head and values picked out for GLG(red), GVG(green) and GHG(blue)');
    plot(t([1 end])',[GLG GLG],'r--');
    plot(t([1 end])',[GVG GVG],'g--');
    plot(t([1 end])',[GHG GHG],'b--');
    
    for i=1:length(HY)
        plot(t(HY(i).J),  h(HY(i).J), 'bo');
        plot(t(HY(i).J(HY(i).GHGJ)),HY(i).GHG,'ro','markerFaceColor','r');
        plot(t(HY(i).K)            ,HY(i).GVG,'go','markerFaceColor','g');
        plot(t(HY(i).J(HY(i).GLGJ)),HY(i).GLG,'bo','markerFaceColor','b');
    end
    plot(t,h,'k');
    datetick('x',12);

    legend('GLG','GVG','GHG','GLG','GVG','GHG','data');
end

%% Show GLG, GVG and GHG for all cross sections
