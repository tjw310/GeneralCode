function estimates = fit180calibs(obj,x1,z1,y1,x2,z2)
%y1=-y1;
%y2=-y2;
st = [-0.1,22.8,-0.1/180*pi,0];
model = @fun;

options = optimset('maxfuneval',2e4,'maxiter',2e4,'tolx',1e-12,'tolfun',1e-12);

estimates = fminsearch(model,st,options);

[~,xplot,zplot,medDist,x,z,rx,rz,ry]= model(estimates);

figure; hold on;
scatter(obj.opticCentre(1),obj.opticCentre(2),'k');
scatter(x1,z1,[],y1,'filled'); scatter(x,z,'rx'); hold off;
colormap(jet);

figure; scatter(xplot,zplot,'rx'); hold on; scatter(x2,z2,'b'); hold off;
title(num2str(medDist));

figure; hold on;
scatter(obj.opticCentre(1),obj.opticCentre(2),'k');
scatter(xplot,zplot,[],ry,'filled'); scatter(rx,rz,'rx'); hold off;
colormap(jet);

function [sse,x2g,z2g,sse2,x,z,rx,rz,ry] = fun(params)
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = params(1);
    phi = params(3);
    gamma = params(4);
    
    %u=-0.1;
    %v=0;
    %w=1;
    
    xa = params(2);
    %ya = params(3);
    ya = -0.009;
    
   % mt = params(6);
    
   % phi =0;
    %gamma = 0;
    
    %xa = 0;
    %ya = -0.006;
    

    y = y1;
    x = (x1-xo).*(1-dmdy*y)+xo;
    z = (z1-zo).*(1-dmdy*y)+zo;
    
    
    u = sin(phi)*cos(gamma);
    v = sin(phi)*sin(gamma);
    w = cos(phi);
    
    f = u*(x-xa)+v*(y-ya)+w*z;
    
    rx = 2*u*f-x+2*xa;
    ry = 2*v*f-y+2*ya;
    rz = 2*w*f-z;
    
    x2g = (rx-xo)./(1-dmdy*ry)+xo;
    z2g = (rz-zo)./(1-dmdy*ry)+zo;
    y2g = ry;
    
    %sse = nansum(sqrt((x2g-x2).^2+(z2g-z2).^2+(y2g-y2).^2));
    %sse2 = nanmedian(sqrt((x2g-x2).^2+(z2g-z2).^2+(y2g-y2).^2))
    
    sse = nansum(sqrt((x2g-x2).^2+(z2g-z2).^2));
    sse2 = nanmedian(sqrt((x2g-x2).^2+(z2g-z2).^2));
    
    scatter(x2,z2); hold on; scatter(x2g,z2g,'rx'); hold off; drawnow;
    %scatter(-x1+xa,z1,[],y1,'filled'); hold on; scatter(x2g,z2g,'rx'); scatter(x2g,z2g,[],y2g,'filled');hold off; drawnow;
end

end