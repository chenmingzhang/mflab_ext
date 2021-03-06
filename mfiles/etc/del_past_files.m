function del_past_files(basename)
% delete the past simulated files 
% now included RCH file
delete([basename,'.BAS'],[basename,'.BCF'],[basename,'.BGT'],...
       [basename,'.DDN'],[basename,'.DIS'],[basename,'.EVT'],...
       [basename,'.HDS'],[basename,'.MF2000.LST'],[basename,'.OC'],...
       [basename,'.PCG'],[basename,'.WEL'],[basename,'.CHD'],...
       [basename,'.WEL'],...
       [basename,'.LPF'],[basename,'.MF2005.LST'],['mf2005.nam'],...
       ['mf2k.nam'],[basename,'.RCH'],[basename,'.DE4'])
