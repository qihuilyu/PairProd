 r = vrrotvec([1,0,0],[0.743145 -0.669131 0.000000]);
 
r(4)/pi*360

r = [0,0,-1,pi/2];
m = vrrotvec2mat(r);
(m*[1,0,0]')'
