load('D:\datatest\PairProd\GBMHY_100m\PairProd\results\GBMHY_100m PairProd Info7 eta0 result.mat')
xPolish = result.xPolish;

xPolish(xPolish<100) = 100;

ds = 25/2*10;
mkdir(['D:\datatest\PairProd\GBM_final\'])
fileID = fopen(['D:\datatest\PairProd\GBM_final\nparticles_GBMHY_100m PairProd Info7 eta0 result.txt'],'w');
fprintf(fileID,'%d \n',ceil(xPolish/ds));
fclose(fileID);


