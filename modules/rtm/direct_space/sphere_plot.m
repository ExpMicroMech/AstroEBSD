screen_int
[x,y,z]=sphere(500);
ints=Cube_Sample(x,y,z,screen_int,0);
figure
surf(x,y,z,rescale(reshape(ints,size(x))),"EdgeColor","none")
colormap gray;
axis equal;
%%
figure
imagesc(reshape(ints,size(x)))

%%
figure
% scatter3(x,y,z,rescale(reshape(ints,size(x)))
surf(x,y,z)