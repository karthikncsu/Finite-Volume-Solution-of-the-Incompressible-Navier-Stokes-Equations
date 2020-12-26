clear all

[x,y] = meshgrid(0:1:20,-3:.5:3);

[m,n] = size(y);
p = 2000;
l = 20;
c = 3;
for i = 1:m
    for j = 1:n
        
       sig11(i,j) = ((0.3*p*y(i,j)* x(i,j))/(l*c))+((p/(4*l*c^3))*y(i,j)*x(i,j)^3)-((p/(2*l*c^3))*x(i,j)*y(i,j)^3);
       %sig22(i,j) = ((-0.5*p*x(i,j))/l)-((3*p/(4*l*c))*y(i,j)*x(i,j))+((p/(4*l*c^3))*x(i,j)*y(i,j)^3);
       %sig12(i,j) = ((p*c)/(40*l))+((3*p/(8*l*c))*x(i,j)^2)-((3*p/(20*l*c))*y(i,j)^2)-((3*p/(8*l*c^3))*x(i,j)^2*y(i,j)^2)+((p/(8*l*c^3))*y(i,j)^4);
    end
end
[d,g] = contourf(x,y,sig11);
clabel(d,g);
xlabel('x')
ylabel('y')
title('sigma11      P = 2000 , l = 20, c = 3')
