
filename = 'C:\Users\admin\Desktop\mcgeo.txt';
A = readmatrix(filename);


density = reshape(A(:,1),[105,125,65]);
figure;imshow3D(density,[])
Ind = reshape(A(:,2),[105,125,65]);
figure;imshow3D(Ind,[])


