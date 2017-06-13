function [ xout,yout,zout,u,v,w ] = rotArbAxis(x,y,z,phi,gamma,t )
%rotation angle t, of point [x;y;z] about arbitary axis with unit vector
%[u;v;w], expressed by angle phi to the z-axis, and angle gamma around.

%phi goes from 0-180. %gamma goes from 0-360

u = cos(phi/180*pi)*sin(gamma/180*pi);
v = sin(phi/180*pi)*sin(gamma/180*pi);
w = cos(phi/180*pi);

t = t/180*pi;

p = [x;y;z;1];

M = [u^2+(1-u^2)*cos(t),u*v*(1-cos(t))-w*sin(t),u*w*(1-cos(t))+v*sin(t),0;...
    u*v*(1-cos(t))+w*sin(t),v^2+(1-v^2)*cos(t),v*w*(1-cos(t))-u*sin(t),0;...
    u*w*(1-cos(t))-v*sin(t),v*w*(1-cos(t))+u*sin(t),w^2+(1-w^2)*cos(t),0;...
    0,0,0,1];

V = M*p;

xout = V(1); yout = V(2); zout = V(3);

end

