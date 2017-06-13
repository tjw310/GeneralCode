function est = fitTrace(xrealArray,zrealArray,beadno,fixed,obj)
t = 360-obj.theta;
dmdy = obj.e(2);
xo = obj.opticCentre(1);
zo = obj.opticCentre(2);

xreal = xrealArray(beadno,:);
zreal = zrealArray(beadno,:);

yGuess = fixed(beadno).zDepth./obj.pxSz.*obj.mtz0;
xGuess = (xreal(1)-xo).*(1-dmdy.*yGuess./obj.mtz0*obj.pxSz)+xo;
zGuess = (zreal(1)-zo).*(1-dmdy*yGuess./obj.mtz0*obj.pxSz)+zo;

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
    dx = 0;
    dy = 0;
    
    for i=1:length(xreal)
    [xp(i),yp(i),zp(i)] = ARtranslate(x,y,z,-dx,-dy);
    [xp(i),yp(i),zp(i)] = rotArbAxis(x,y,z,phi,gamma,t(i));
    [xp(i),yp(i),zp(i)] = ARtranslate(xp(i),yp(i),zp(i),dx,dy);
    [xp(i),yp(i),zp(i)] = magScale(xp(i),yp(i),zp(i),dmdy.*obj.pxSz/obj.mtz0,xo,zo);
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