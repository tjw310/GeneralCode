function z = testFanv2(obj,projNo)
u = obj.x-obj.opticCentre(1)/obj.SubSampFc;

beta = linspace(0,2*pi-2*pi/obj.nProj,obj.nProj);
beta= beta(projNo);

sqrd = obj.opticCentre(1)/obj.SubSampFc;
OaDl = obj.nPx/2-abs(obj.AoRcentreX/obj.SubSampFc);
OsOa = (obj.AoRcentreX-obj.opticCentre(1))/obj.SubSampFc;
OsDl = OsOa+OaDl;
D = obj.D;

prefact = (obj.D-u./obj.D.*OsOa)./sqrt(obj.D.^2+u.^2);

proj = squeeze(obj.projections(:,size(obj.projections,1)/2,projNo)).*prefact.';

gamma = atan2(u,D);
theta = atan2(OsDl,D);
rho = gamma - atan2(OsOa,D);
ep = atan2(OsDl,D)-atan2(OsOa,D);

weight = 0.5*(sin(pi*rho./(-2*ep))+1);
weight(rho<-ep)=1;

p = proj.*weight.';

filt_len = max(64,2^nextpow2(2*obj.nPx-1)); 
lBefore = round((filt_len-obj.nPx)/2-obj.opticCentre(1)/obj.SubSampFc);
lAfter = filt_len-obj.nPx-lBefore;

%lBefore = round((filt_len-obj.nPx)/2); lAfter = round((filt_len-obj.nPx)/2);

p = padarray(p,[lBefore,0],'pre','replicate');
p = padarray(p,[lAfter,0],'post','replicate');

ramp = abs(linspace(-1,1,size(p,1)).');

f = fft(p).*fftshift(ramp);

g = ifft(f);
g = g(lBefore+1:lBefore+obj.nPx);

[xx,yy]= meshgrid(obj.x,obj.y);
pos = (obj.AoRcentreX+obj.AoRmotion(projNo))/obj.SubSampFc;
l = obj.D.*(xx.*cos(beta)+yy.*sin(beta)-sqrd+pos)./(xx.*sin(beta)-yy.*cos(beta)+obj.D)+sqrd;

g = repmat(g.',size(l,1),1);

z = interp2(xx,yy,g,l,yy);

%z(isnan(z)) = 0;

z = z.*obj.D^2./(xx*sin(beta)-yy.*cos(beta)+obj.D).^2/2;





end