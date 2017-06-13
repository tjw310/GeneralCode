function testTraceDif(phi,gamma,motx,moty,images)

dmdy = images.e(2)*images.pxSz/images.mtz0;
xo = images.opticCentre(1);
zo = images.opticCentre(2);

t = 360-images.theta;

clear xi yi zi xp yp zp x y z

testx = [xo+200,xo-200,xo+200,xo-200,xo+200,xo-200,xo+200,xo-200,xo+200,xo-200,xo+200,xo-200];
testy = [200,200,200,200,-200,-200,-200,-200,0,0,0,0];
testz = [zo+200,zo+200,zo-200,zo-200,zo+200,zo+200,zo-200,zo-200,zo+200,zo+200,zo-200,zo-200];

%h = figure;
h2 = figure;

for i=1:length(testx)
for k=1:length(t)
    angle = t(k);
    
    [xi(k),yi(k),zi(k),u,v,w] = rotArbAxis(testx(i),testy(i),testz(i),0,0,angle); 
    [xi(k),yi(k),zi(k)] = magScale(xi(k),yi(k),zi(k),dmdy,xo,zo);
    
    
    [x,y,z] = rotArbAxis(testx(i),testy(i),testz(i),phi,gamma,angle);
    [xp(k),yp(k),zp(k)] = ARtranslate(x,y,z,motx(k),moty(k));
    [xp(k),yp(k),zp(k)] = magScale(xp(k),yp(k),zp(k),dmdy,xo,zo);
    
    xShift(k) = xi(k)-xp(k);
    zShift(k) = zi(k) - zp(k);
end

%figure(h);
%plot3(xi,yi,zi); 
%hold on; scatter3(xp,yp,zp,[],1:400,'filled'); hold off;

figure(h2);
subplot(3,4,i)
plot(xShift); hold on; plot(zShift); hold off; %legend('xshift','zshift');
title(num2str([testx(i)-xo,testy(i),testz(i)-zo],1));

end
drawnow;
end