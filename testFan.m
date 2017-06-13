function z = testFan(obj,projNo)

u = obj.x-obj.opticCentre(1)/obj.SubSampFc;

beta = linspace(0,2*pi-2*pi/obj.nProj,obj.nProj);
beta= beta(projNo);

pos = (obj.AoRcentreX+obj.AoRmotion(projNo))/obj.SubSampFc;
sqrd = obj.opticCentre(1)/obj.SubSampFc;

prefact = (obj.D+u./obj.D.*(sqrd-pos))./sqrt(obj.D.^2+u.^2);

proj = squeeze(obj.projections(:,size(obj.projections,1)/2,projNo)).*prefact.';

p = proj;

filt_len = max(64,2^nextpow2(2*obj.nPx-1)); 
lBefore = round((filt_len-obj.nPx)/2+obj.opticCentre(1)/obj.SubSampFc);
lAfter = filt_len-obj.nPx-lBefore;

p = padarray(p,[lBefore,0],'pre','replicate');
p = padarray(p,[lAfter,0],'post','replicate');

ramp = abs(linspace(-1,1,size(p,1)).');

f = fft(p).*fftshift(ramp);

g = ifft(f);
g = g(lBefore+1:lBefore+obj.nPx);

[xx,yy]= meshgrid(obj.x,obj.y);

l = obj.D.*(xx.*cos(beta)+yy.*sin(beta)-sqrd+pos)./(xx.*sin(beta)-yy.*cos(beta)+obj.D)+sqrd;

g = repmat(g.',size(l,1),1);

z = interp2(xx,yy,g,l,yy);

%z(isnan(z)) = 0;

z = z.*obj.D^2./(xx*sin(beta)-yy.*cos(beta)+obj.D).^2/2;


end