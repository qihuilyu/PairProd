nx = 60;
ny = 5;
nz = 20;

mask0 = zeros(nx,ny,nz);
mask0(56:57,2:3,1) = 1;
figure;imshow3D(mask0,[]);

airHU = -1000;
waterHU = 0;

phantom = ones(nx,ny,nz)*waterHU;
phantom(mask0==1) = 3000;

figure;imshow3D(phantom,[-1000,3000])

BODY = ones(nx,ny,nz);
BODY(:,:,10:20) = 0;
BODY(1:40,:,10:20) = 0;
phantom(BODY==0) = airHU;
PTV = mask0;
figure;imshow3D([BODY PTV],[])

save('D:\datatest\PairProd\phantom_dimtest\phantom_dimtest.mat','phantom','BODY','PTV')









