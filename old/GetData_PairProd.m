%% Dicom
DicomPath = '/media/raid1/qlyu/PairProd/datatest/phantom/data/water_phantom';
baseFileNames = dir([DicomPath '/*.dcm']);
[sortedFile,index] = sort_nat({baseFileNames.name});

img = [];
for ii = 1:numel(sortedFile)
    img(:,:,ii) = dicomread(fullfile(DicomPath,sortedFile{ii}));
end
figure;imshow3D(img,[0,2000])

%%
LIsize = [10 10 10];
filename = '/media/raid1/qlyu/VHEE/code_QL/montecarlo/usercode/code/build/NofPositronAnni3D.bin';
fid = fopen(filename);
FPx=fread(fid,LIsize(1)*LIsize(2)*LIsize(3),'double');
fclose(fid);

sum(FPx)
FPx = reshape(FPx,LIsize);
figure;imshow3D(FPx,[])



filename = '/media/raid1/qlyu/VHEE/code_QL/montecarlo/usercode/code/build/ids.csv';
M = csvread(filename,1,0);
DetectorID = M(:,1);
Time = M(:,2);
Energy = M(:,3);
EventID = M(:,4);
ThreadID = M(:,5);

Ind_511 = find(abs(Energy-0.511)<0.0001);
NumAnniPhotons = (length(Ind_511)-length(unique(EventID(Ind_511))))*2;

numphotons = size(M,1);
NumAnniPhotons/numphotons

EnergyResolution = 0.1;
Ind_accept = find(abs(Energy-0.511)<0.511*EnergyResolution);
NumAcceptedPhotons = length(Ind_accept);

Ind_10MeV = find(abs(Energy-10)<0.01);
NumPrimary = length(Ind_10MeV);


doserate = 0.1/60; % (0.1Gy/min)
numevent = 1e+06;
time = max(dose)*numevent/doserate;
eventrate = time/numevent*1e+09; % ns

deltatime = normrnd(eventrate,eventrate/5,numevent,1);
% figure;hist(deltatime)
cumsumtime = cumsum(deltatime);
CorrectedTime = Time + cumsumtime(EventID+1);
[CorrectedTime(Ind_511) M(Ind_511,:)]

CoincidenceTime = 0.02;  % ns
[sortedtime, sortInd] = sort(CorrectedTime);
sortInd_coin = find(diff(sortedtime)<CoincidenceTime);
Ind_coin1 = sortInd(sortInd_coin);
Ind_coin2 = sortInd(sortInd_coin+1);

% test = [M(Ind_coin1,:) M(Ind_coin2,:)];
% test2 = [CorrectedTime(Ind_coin1) CorrectedTime(Ind_coin2)]
[sortedtime(sortInd_coin) sortedtime(sortInd_coin+1)]

[Ind_coin1_511buff, iInd_coin1, iInd_511_1] = intersect(Ind_coin1, Ind_511);
[Ind_coin2_511buff, iInd_coin2, iInd_511_2] = intersect(Ind_coin2, Ind_511);

iInd_coin = intersect(iInd_coin1,iInd_coin2);
Ind_coin1_511 = Ind_coin1(iInd_coin);
Ind_coin2_511 = Ind_coin2(iInd_coin);

test = [M(Ind_coin1_511,:) M(Ind_coin2_511,:)];


[Ind_coin1_acceptbuff, iInd_coin1, iInd_accept_1] = intersect(Ind_coin1, Ind_accept);
[Ind_coin2_acceptbuff, iInd_coin2, iInd_accept_2] = intersect(Ind_coin2, Ind_accept);

iInd_coin = intersect(iInd_coin1,iInd_coin2);
Ind_coin1_accept = Ind_coin1(iInd_coin);
Ind_coin2_accept = Ind_coin2(iInd_coin);

test = [M(Ind_coin1_accept,:) M(Ind_coin2_accept,:)];


length(Ind_coin1_511)
length(Ind_coin2_511)*2
length(Ind_coin1)
length(Ind_coin1)*2
length(Ind_coin1_accept)
length(Ind_coin1_accept)*2
