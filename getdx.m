function [dx,dz] = getdx( xrealArray,zrealArray,obj,fixed )
% gets dx and dz based on first xreal,zreal.

t = 360-obj.theta;
dmdy = obj.e(2);
xo = obj.opticCentre(1);
zo = obj.opticCentre(2);

for i=1:size(xrealArray,1)
    yGuess(i) = fixed(i).zDepth./obj.pxSz.*obj.mtz0;
    xGuess(i) = (xrealArray(i,1)-xo).*(1-dmdy.*yGuess(i)./obj.mtz0*obj.pxSz)+xo;
    zGuess(i) = (zrealArray(i,1)-zo).*(1-dmdy*yGuess(i)./obj.mtz0*obj.pxSz)+zo;


for k=1:length(t);
    [x,y,z] = rotArbAxis(xGuess(i),yGuess(i),zGuess(i),0,0,t(k));
    [xout(k),yout(k),zout(k)] = magScale( x,y,z,dmdy*obj.pxSz/obj.mtz0,xo,zo );
end

dx(i,:) = xrealArray(i,:)-xout;
dz(i,:) = zrealArray(i,:)-zout;

figure; plot(dx(i,:)); hold on; plot(dz(i,:)); hold off; drawnow;
disp('done');
%figure;
%scatter(xrealArray(i,:),zrealArray(i,:)); hold on; scatter(xout,zout,'rx'); hold off;
%drawnow;

end



