%% Example see USGS modpath Version 6 (2012)
%% SIMULATION3: backward endpoint simulation (see manual p36ff 43)
% TO 130221

path(path,fileparts(pwd))

%% make sure path is set to parent directory
d=dir('..');
if ~strmatchi('mf_adaptALL.m',{d.name})
    error('%s: Can''t find file %s in the parent directory,\n%s\n',...
        mfilename,'mf_adaptALL.m',fileparts(pwd));
else        
    mf_adaptALL;
end

%% Modpath info pertaining to this simulation

%%
% Specify the number and placement of the starting points with each cell.
% in this case cells are in top of all rows, column 5 layer 1.
zone=[15 20 8 12 2 4]; % specify the cell range and don't use ZONE

zoneVals = {zone,'name',basename,'placement',2,'LineSpec','bo'}; 

%%
% Generate the mpath_particleGroupObj from which MODPATH can generate the
% particle starting locations

pGrp = mpath_particleGroupObj(gr,zoneVals);

%% Alternative
%pGrp = mpath_particleGroupObj(gr,ZONE,zone,'name',basename,'IFace',[1 2 3 4 6],'placement',10,'LineSpec','bo');

%% Particles
% The definition above will allow to generate input for MODPATH from which
% MODPATH can generate the required starting points. So mfLab does not need
% to do that. To allow plotting the particles within mfLab, they can also
% be generated by mfLab using the method getParticles. This methods addes
% the particles to the mpath_particleGroupObj
pGrp   = pGrp.getParticles(gr);

%% Show particles in 3D

figure; hold on; view(3); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

gr.plotMesh('faceAlpha',0.15);

pGrp.plot(); title('Particles starting points');

%% You can turn the graphic by hand to better view the particles

save underneath zoneVals
