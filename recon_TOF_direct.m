function img_direct = recon_TOF_direct(CorrectedTime, detectorIds, Ind_coin_accept, cilist, reconparams)

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
n_p1_p2 = p1_p2./sqrt(sum(p1_p2.^2,2));
p0 = -n_p1_p2.*deltar + (p1+p2)/2;

%% Show image
p0(sqrt(sum(p0.^2,2))>distrange) = NaN;
cilist(isnan(p0(:,1))) = [];
p0(isnan(p0(:,1)),:)=[];
xp0 = p0(:,1);
yp0 = p0(:,2);

img_direct = zeros(imgsize(1),imgsize(2));
xImage = ((1:imgsize(1))-(1+imgsize(1))/2)*imgres;
yImage = ((1:imgsize(2))-(1+imgsize(2))/2)*imgres;

for ii = 1:length(xp0)
    ixp0 = xp0(ii);
    iyp0 = yp0(ii);
    
    if(ixp0<min(xImage)||ixp0>max(xImage)||iyp0<min(yImage)||iyp0>max(yImage))
         continue
    end
    
    [~,xind] = min(abs(ixp0-xImage));
    [~,yind] = min(abs(iyp0-yImage));

    img_direct(xind,yind) = img_direct(xind,yind) + 1/cilist(ii);
end





