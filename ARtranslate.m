function [ xout,yout,zout ] = ARtranslate( x,y,z,arx,ary )
%translates point x,y,z using vector [arx;ary;0];

p = [x;y;z;1];

T = [1,0,0,arx;0,1,0,ary;0,0,1,0;0,0,0,1];

V = T*p;

xout = V(1); yout = V(2); zout = V(3);

%scatter(x,z); hold on; scatter(xout,zout); hold off;
%drawnow; pause(1);


end

