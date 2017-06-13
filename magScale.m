function [ xout,yout,zout ] = magScale( x,y,z,dmdy,xo,zo )
%scales x,z coordinates based on their y position and the dmdy factor,
%optic centre coords xo,zo

xout = (x-xo)/(1-dmdy*y)+xo;
zout = (z-zo)/(1-dmdy*y)+zo;
yout = y;


end

