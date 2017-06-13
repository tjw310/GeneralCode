function estimates = fit180projections(obj,fixed)
[xinit,zinit] = FastPeakFind(obj.projections(:,:,1).');
[x2,z2] = FastPeakFind(obj.projections(:,:,201).');

k=1;
for i=1:length(xinit)
    for j=1:length(fixed)
        r(j) = sqrt((fixed(j).centreX-xinit(i))^2+(fixed(j).centreY-zinit(i))^2);
    end
    [mn,lc] = min(r);
    if mn<75
        y1(k) = fixed(lc).zDepth;
        x1(k) = xinit(i)-2560/2;
        z1(k) = zinit(i)-2160/2;
        k=k+1;
    end
end

st = [-0.1,20,0,0,0];
model = @fun;

options = optimset('maxfuneval',2e4,'maxiter',2e4,'tolx',1e-14,'tolfun',1e-14);

estimates = fminsearch(model,st,options);

[~,xplot,zplot] = model(estimates);

figure;
imagesc(obj.projections(:,:,201).'); hold on;
scatter(x2,z2,'w'); scatter(xplot+2560/2,zplot+2160/2,'rx'); hold off; drawnow;

function [sse,x2gcopy,z2gcopy] = fun(params)
    x2copy = x2-2560/2;
    z2copy = z2-2160/2;
    sse = 0;
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = params(1);
    phi = params(4)/180*pi;
    gamma = params(5)/180*pi;
    
    %phi = 0;
    %gamma = 0;
    
   xa = params(2);
   ya = params(3);
    
    theta =pi;    

    y = y1;
    x = (x1-xo).*(1-dmdy*y)+xo;
    z = (z1-zo).*(1-dmdy*y)+zo;
    
    
    u = sin(phi)*cos(gamma);
    v = sin(phi)*sin(gamma);
    w = cos(phi);
    
    f = (u*(x-xa)+v*(y-ya)+z)*(1-cos(theta));
    
    rx = u*f+(x-xa)*cos(theta)+(-w*(y-ya)+v*z)*sin(theta)+xa;
    ry = v*f+(y-ya)*cos(theta)+(w*(x-xa)-u*z)*sin(theta)+ya;
    rz = w*f+z*cos(theta)+(-v*(x-xa)+u*(y-ya))*sin(theta);
    
    x2g = (rx-xo)./(1-dmdy*ry)+xo;
    z2g = (rz-zo)./(1-dmdy*ry)+zo;
    y2g = ry;
    
    x2gcopy = x2g;
    z2gcopy = z2g;
    
   % subplot(1,2,1); 
    %scatter(x2copy,z2copy); hold on;
    %subplot(1,2,2);
    %scatter(x2g,z2g,'rx'); hold off; drawnow;
    
    l=1;
    while ~isempty(x2copy)
        r2 = sqrt((x2g-x2copy(1)).^2+(z2g-z2copy(1)).^2);
        [mn2,lc2] = min(r2);
        allsse(l) = mn2(1);
        x2copy = x2copy(2:end);
        z2copy = z2copy(2:end);
        l=l+1;
        x2g = x2g((1:length(x2g))~=lc2);
        z2g = z2g((1:length(z2g))~=lc2);
    end
    %sse
    sse = nansum(allsse);
    nanmedian(allsse)
    
    
end

end