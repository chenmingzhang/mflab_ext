% Testing Piez....

time = [Hds.time];  % PiezM1L2Q65(1).UserData.time;

Np     = numel(PiezM1L2Q65);

ddM1L2 = NaN(Np,numel(time));
ddM2L2 = NaN(Np,numel(time));

for ip = 1:Np
   
    ddM1L2(ip,:) = PiezM1L2Q65(ip).UserData.hd - PiezM1L2Q00(ip).UserData.hd;
    ddM2L2(ip,:) = PiezM2L2Q65(ip).UserData.hd - PiezM2L2Q00(ip).UserData.hd;
end

figure('name','dd vs time','pos',screenPos(1,0.6));
hold on;
for ip=1:Np
    plot(time,ddM1L2(ip,:),[mf_color(ip),'-'],'lineWidth',2);
    plot(time,ddM2L2(ip,:),[mf_color(ip),'.'],'lineWidth',1);
end