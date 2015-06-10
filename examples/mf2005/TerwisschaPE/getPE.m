%% Get KNMI data for station Leeuwarden

if ~exist('KNMI.mat','file');
    readKNMI('KNMI_20131209.txt'); % saves KNMI.mat holding knmiData and knmiMeta    
end

load KNMI  % for Leeuwarden

iP = strmatchi('RH',knmiMeta(:,1),'exact');
iE = strmatchi('EV24',knmiMeta(:,1),'exact');
iT = strmatchi('YYYY',knmiMeta(:,1));

% Convert time to datenums
T    = knmiData(:,iT);
yy   = floor(T/10000);
mmdd = (T-10000*yy);
mm   = floor(mmdd/100);
dd   = mmdd-mm*100;
TPE  = [datenum(yy,mm,dd), knmiData(:,[iP iE])/10];

%% change neg values to zero and remove periods with NaN values
TPE(TPE<0)=0;
TPE = TPE(~isnan(TPE(:,end)),:);

%% Show what we got
figure; hold on;

title(sprintf('Precip and Makkink for STN %d [%s]',knmiData(1,1), knmiMeta{1,2}));
plot(TPE(:,1),TPE(:,2),'b');
plot(TPE(:,1),TPE(:,3),'g');
xlabel('time');
ylabel('mm/d');
datetick();

%% Select a period of 10 years

TPEE = TPE(TPE(:,1)<TPE(1,1)+10*365.24);

DV = [datevec(TPE(:,1)) TPE(:,2:3)/1000];
DV = DV(:,[1:3,end-1:end]);


fprintf('%d\t%d\t%d\t%g\t%g\n',DV');



