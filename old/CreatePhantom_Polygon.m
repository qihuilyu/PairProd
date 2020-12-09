nx = 141;
ny = 153;
nz = 100;

x = [13 125 136 121 54];
y = [10 10  105 136 114];
bw = poly2mask(x,y,nx,ny);
figure;imshow(bw)

airHU = -1000;
waterHU = 0;
HUs = [-1000,-400,-150,100,300,2000,4927,6000];

[x,y,z] = ndgrid(1:nx,1:ny,1:nz);

phantom = ones(nx,ny,nz)*airHU;
mask0 = repmat(bw,[1,1,nz]);
phantom(mask0) = waterHU;

inthu = [-5000.0, -1000.0, -400, -150, 100, 300, 2000, 4927, 66000];
intdens = [0.0, 0.01, 0.602, 0.924, 1.075, 1.145, 1.856, 3.379, 7.8];
huq = -5000:66000;
densq = interp1(inthu,intdens,huq);
figure;plot(huq,densq)

matdens = [0.0,1.654, 7];
for ii = 1:numel(matdens)
    [~,ind] = min((matdens(ii)-densq).^2);
    mathu(ii) = huq(ind);
end

r = 30; ninserts = numel(matdens); ri = 8;
for ii = 1:ninserts
    theta = ii/ninserts*2*pi + pi/6;
    rx = (1+nx)/2 + r*cos(theta)-10;
    ry = (1+ny)/2 + r*sin(theta)+10;
    mask = ((x-rx).^2+(y-ry).^2 <ri.^2 & mask0==1);
    phantom(mask) = mathu(ii);
end

figure;imshow3D(phantom,[-1000,3000])

BODY = mask0;
PTV = mask0;
PTV(:,:,[1:48 53:100]) = 0;
figure;imshow3D([BODY PTV],[])

save('D:\datatest\PairProd\phantom_polygon\phantom_polygon.mat','phantom','BODY','PTV')









