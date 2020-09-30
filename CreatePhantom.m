
nx = 301;
ny = 301;
nz = 301;
airHU = -1000;
waterHU = 0;
HUs = [-1000,-400,-150,100,300,2000,4927,6000];

[x,y,z] = ndgrid(1:nx,1:ny,1:nz);

phantom = ones(nx,ny,nz)*airHU;
mask0 = ((x-(1+nx)/2).^2/(120.^2)+(y-(1+ny)/2).^2/(100.^2)<1 & z<280 & z>20);
phantom(mask0) = waterHU;

inthu = [-5000.0, -1000.0, -400, -150, 100, 300, 2000, 4927, 66000];
intdens = [0.0, 0.01, 0.602, 0.924, 1.075, 1.145, 1.856, 3.379, 7.8];
huq = -5000:66000;
densq = interp1(inthu,intdens,huq);
figure;plot(huq,densq)

matdens = [0.0,0.207,0.481,0.919,0.979,1.004,1.109,1.113,1.496,1.654];
for ii = 1:numel(matdens)
    [~,ind] = min((matdens(ii)-densq).^2);
    mathu(ii) = huq(ind);
end

r = 70; ninserts = numel(matdens); ri = 10;
for ii = 1:ninserts
    theta = ii/ninserts*2*pi;
    rx = (1+nx)/2 + r*cos(theta);
    ry = (1+ny)/2 + r*sin(theta);
    mask = ((x-rx).^2+(y-ry).^2 <ri.^2 & mask0==1);
    phantom(mask) = mathu(ii);
end

figure;imshow3D(phantom,[-1000,1600])

save('D:\datatest\PairProd\phantom\phantom.mat','phantom')








