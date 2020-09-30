function [sino, inter, Sino3Dbuff] = rebin_PET_TOF(CorrectedTime, detectorIds, Ind_coin_accept, reconparams)

nb_cryst = reconparams.nb_cryst;
R1 = reconparams.R1;
distrange = reconparams.distrange;
imgres = reconparams.imgres;
imgsize = reconparams.imgsize;

Ind_coin1_accept = Ind_coin_accept(:,1);
Ind_coin2_accept = Ind_coin_accept(:,2);
detid_pair = [detectorIds(Ind_coin1_accept) detectorIds(Ind_coin2_accept)];

deltat = CorrectedTime(Ind_coin1_accept) - CorrectedTime(Ind_coin2_accept);
clight = 300; % c = 300mm/ns
deltar = deltat*clight/2; 

theta = detid_pair/nb_cryst*2*pi;
x_ = R1*sin(theta);
y_ = R1*cos(theta);

p1 = [x_(:,1), y_(:,1)];
p2 = [x_(:,2), y_(:,2)];
p1_p2 = p1 - p2;
pC = (p1+p2)/2;
n_p1_p2 = p1_p2./sqrt(sum(p1_p2.^2,2));
p0 = -n_p1_p2.*deltar + (p1+p2)/2;

%% Show image
p0(sqrt(sum(p0.^2,2))>distrange) = NaN;

maskp0 = isnan(p0(:,1));
pC(maskp0,:) = [];
p0(maskp0,:)=[];
p1_p2(maskp0,:)=[];

dist = round(sqrt(sum(pC.^2,2)),4);
xdy = p1_p2(:,1)./p1_p2(:,2);
ang = mod(round(atan(xdy),4),3.1416);

[uniang, ~, indang] = unique(ang);

xC = pC(:,1);
yC = pC(:,2);
dist(yC<0) = -dist(yC<0);
dist(yC==0) = dist(yC==0).*sign(xC(yC==0));
[unidist, ~, inddist] = unique(dist);

maxunidist = max(unidist);
minunidist = min(unidist);
inter = (maxunidist - minunidist)/(numel(unidist)-1);
newunidist = (minunidist:inter:maxunidist)';

p0_pC = p0-pC;
xp0pC = p0(:,1);
yp0pC = p0(:,2);
distp0 = round(sqrt(sum((p0_pC).^2,2)),4);
distp0(xp0pC<0) = -distp0(xp0pC<0);
distp0(xp0pC==0) = distp0(xp0pC==0).*sign(yp0pC(xp0pC==0));

nd = numel(unidist);
na = numel(uniang);
Sino3Dbuff = zeros(nd,nd,na);

TR = 0.02;
sigma = clight*TR/2;
gaussFilter = exp(-(unidist'-distp0) .^ 2 / (2 * sigma ^ 2));
dx = unidist'-distp0;
gaussFilter(dx<-3*sigma|dx>3*sigma) = 0;
gaussFilter = gaussFilter./sum (gaussFilter,2); % normalize
% figure;plot(gaussFilter(12355,:))

for ii = 1:length(distp0)
    id = inddist(ii);
    ia = indang(ii);
    ikernel = gaussFilter(ii,:);
    if(~any(isnan(ikernel)))
        Sino3Dbuff(id,:,ia) = Sino3Dbuff(id,:,ia) + ikernel;
    end
end

sumbuff1 = sum(sum(Sino3Dbuff(1:2:end,:,1)));
sumbuff2 = sum(sum(Sino3Dbuff(2:2:end,:,1)));

sino = 0*Sino3Dbuff;
for ii = 1:na
    for jj = 1:nd
        Ind = find(Sino3Dbuff(:,jj,ii)>0);
        if(numel(Ind)>1)
            vq1 = interp1(unidist(Ind),Sino3Dbuff(Ind,jj,ii),newunidist,'linear');
        elseif(numel(Ind)==1)
            [~,minInd] = min(abs(unidist(Ind)-newunidist));
            vq1 = zeros(nd,1);
            vq1(minInd) = unidist(Ind);
        else
            vq1 = zeros(nd,1);
        end
        % figure(100);plot(unidist(Ind),Sino3Dbuff(Ind,jj,ii),'o',newunidist,vq1,':.'); hold on;
        sino(:,jj,ii) = vq1;
    end
end
sino(isnan(sino)) = 0;
sino(sino<0) = 0;
% figure;imshow(Sino,[])
