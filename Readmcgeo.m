
filename = '/media/raid1/qlyu/PairProd/dbdata/mcgeometry/5f98ac6da39933656826827a/mcgeo.txt';
A = readmatrix(filename);


density = reshape(A(:,1),[105,125,65]);
figure;imshow3D(density,[])

Ind = reshape(A(:,2),[105,125,65]);
figure;imshow3D(Ind,[])


