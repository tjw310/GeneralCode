function est = fitTrace3(x,z,obj)
t = 360-obj.theta;
dmdy = obj.dmdy*obj.pxSz/obj.mtz0;
xo = obj.opticCentre(1);
zo = obj.opticCentre(2);

xreal = x;
zreal = z;

est = fitSinusoid(x); [estz,fit] = fitSinusoid(z);
yGuess = est(1).*sin(est(2)-pi/2)
if or(isnan(x(1)),isnan(z(1)))
xGuess = est(1).*sin(est(2))+est(3);
zGuess = estz(1).*sin(estz(2))+estz(3);
else
    xGuess = (xreal(1)-xo)*(1-dmdy*yGuess)+xo;
    zGuess = (zreal(1)-zo)*(1-dmdy*yGuess)+zo;
end

st = [xGuess,yGuess,zGuess];
model = @fun;
options = optimset('maxfunevals',2e3*length(st),'maxiter',2e3*length(st),'tolfun',1e-8,'tolx',1e-8);
est = fminsearch(model,st,options);

%fun([0,0]);

function sse = fun(params)
     x = params(1);
     y = params(2);
     z = params(3);
    %x = xGuess;
    %y = yGuess;
    %z = zGuess;
    
    %phi = params(4);
    %gamma = params(5);
    %dx = -5.8;
    %dy = 3;
    phi = 0;
    gamma = 0;
    dx = -45;
    dy = 0;
    
    for i=1:length(xreal)
    [xp(i),yp(i),zp(i)] = ARtranslate(x,y,z,-dx,-dy);
    [xp(i),yp(i),zp(i)] = rotArbAxis(xp(i),yp(i),zp(i),phi,gamma,t(i));
    [xp(i),yp(i),zp(i)] = ARtranslate(xp(i),yp(i),zp(i),dx,dy);
    [xp(i),yp(i),zp(i)] = magScale(xp(i),yp(i),zp(i),dmdy,xo,zo);
    end
    
    xp(or(xreal<=-1018,xreal>=1018)) = NaN;
    zp(or(zreal<=-1018,zreal>=1018)) = NaN;
    
    xShift = xp-xreal;
    zShift = zp-zreal;
    
    %plot(xShift); hold on; plot(zShift); hold off; drawnow;
    
    scatter(xreal,zreal); 
    hold on; 
    scatter(xp,zp); 
    hold off; 
    drawnow;
    
    sse = nanmedian(sqrt((xreal-xp).^2+(zreal-zp).^2));
    %sse = nansum(abs(zp-zreal));
end

end