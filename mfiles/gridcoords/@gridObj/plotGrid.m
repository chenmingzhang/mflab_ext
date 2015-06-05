function h = plotGrid(o,varargin)
%% h = o.plotGrid(ax,clr,well,z,figname,figcoords)
% PLOTGRID: Plots the grid lines in color clr given the coordinates xGr yGr
%    and possibly well locations in blue if well is a struct
%    whose elements have fields x and y
%
% USAGE:
%   griDObj.plotgrid(plotOptions)
%   gridObj.plotGrid('color','k','edgealpha','0.2','lineStyle','--');
%
% SEE ALSO: well.plot3D gridObj.googleMap plotgrid
% TO 091201 100115 121016


default = {'edgecolor','k','edgealpha',0.2,'facecolor','none','visible','on'};

[ax, varargin] = getNext(varargin,'axis',[]);
[ax, varargin] = getProp(varargin,'axis',ax);

% default axis properties
[xlim, varargin] = getProp(varargin,'xlim',o.xGr([1 end]));
[ylim, varargin] = getProp(varargin,'ylim',o.yGr([end 1]));
[   ~, varargin] = getProp(varargin,'zlim',o.zGr([end 1]));

[lineSpec, varargin] = getNext(varargin,'char','c');
[IL,c,~,LS] = isLineSpec(lineSpec);

[well, varargin] = getNext(varargin,'wellObj',[]);
[well, varargin] = getNext(varargin,'MNW1Obj',well);
[well, varargin] = getNext(varargin,'MNW2Obj',well);

[figName,varargin] = getNext(varargin,'char','');
[figName,varargin] = getProp(varargin,'fig',figName);
[figPos ,varargin] = getProp(varargin,'figPos',screenPos(0.75));

if IL
    if ~isempty(c)
        I = strmatchi('edgeColor',default);
        default{I+1}=c;
    end
    if ~isempty(LS)
        default = [default, 'lineStyle', LS];
    end
end

[iz, varargin] = getProp(varargin,'iz',[]);
[iz, varargin] = getNext(varargin,'double',iz);

if isempty(iz)
    [z,varargin] = getProp(varargin,'z',[]);
    [z,   ~    ] = getNext(varargin,'double',z);
    if isempty(z)
        iz=1;
    else
        iz = hit(o.zGr,z);
        if isnan(iz)
            error('%s: z must be between %g and %g',mfilename,o.zGr([1 end]));
        end
    end
else
    if iz<1 || iz>o.Nlay
        error('%s: iz must be between 1 and %d',mfilename,o.Nlay);
    end
end

if isempty(figName)
    set(gcf,'position',figPos);
else
    figure('name',figName,'position',figPos);
end

axProps = {'nextplot','add','xgrid','on','ygrid','on','zgrid','on','xlim',xlim,'ylim',ylim};
if isempty(ax)
    ax = gca;
end

set(ax,axProps{:})

if ~isempty(well)
    well.plotXY(ax);
end

h = surf(o.xGr,o.yGr,o.zGr(iz)*ones(o.Ny+1,o.Nx+1),default{:},'parent',ax);

