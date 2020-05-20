function [x,y,z] = GenGrid(grid_N,wr,wi)
x = linspace(wr(1),wr(2),grid_N);
y = linspace(wi(1),wi(2),grid_N);
[x_,y_] = meshgrid(x,y);
z = x_ + 1i*y_;
end