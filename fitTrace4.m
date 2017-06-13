function [est,motorMotion] = fitTrace4(x,z,obj)
t = obj.theta;
dmdy = obj.dmdy;
xo = obj.opticCentre(1);
zo = obj.opticCentre(2);

xreal = x;
zreal = z;

[est,fit] = fitSinusoid(x);
yGuess = ((max(x)-min(x))/2)*sin(t/180*pi+est(2)+pi/2);
xGuess = ((max(x)-min(x))/2)*sin(t/180*pi+est(2)-pi)+est(3);
zGuess = linspace((max(z)+min(z))/2,(max(z)+min(z))/2,obj.nProj);

xS = (xGuess-xo)./(1-dmdy*obj.pxSz/obj.mtz0*yGuess)+xo;
zS = (zGuess-zo)./(1-dmdy*obj.pxSz/obj.mtz0*yGuess)+zo;

figure; plot(x); hold on; plot(fit); hold off;
figure; scatter(x,z,[],1:obj.nProj); hold on; scatter(xS,zS,[],1:obj.nProj); hold off; 
title(num2str(nanmedian(sqrt((xreal-xS).^2+(zreal-zS).^2)))); drawnow;


figure;

st = [(max(x)-min(x))/2,est(2),est(3),(max(z)+min(z))/2,zo]
model = @fun;
options = optimset('maxfunevals',2e3*length(st),'maxiter',2e3*length(st),'tolfun',1e-8,'tolx',1e-8);
est = fminsearch(model,st,options);
[~,motorMotion] = model(est);

figure; plot(motorMotion);

%fun([0,0]);

function [sse,motorMotion] = fun(params)
    A = params(1);
    phase = params(2);
    dm = obj.dmdy;
    phi = 0;
    gamma = 0;
    dx = params(3);
    dy = 0;
    zo = params(5);
    zval = params(4);
    
    y2 = A*sin(t/180*pi+phase+pi/2);
    x2 = A*sin(t/180*pi+phase-pi)+dx;
    z2 = linspace(zval,zval,obj.nProj);
    
    xp = (x2-xo)./(1-dm*obj.pxSz/obj.mtz0*y2)+xo;
    zp = (z2-zo)./(1-dm*obj.pxSz/obj.mtz0*y2)+zo;
    
    
    %xp(or(xreal<=-1018,xreal>=1018)) = NaN;
    %zp(or(zreal<=-1018,zreal>=1018)) = NaN;
    
    xShift = xp-xreal;
    zShift = zp-zreal;
    
    motorMotion = (xreal-xo).*(1-dm*obj.pxSz/obj.mtz0*y2)+xo-x2;
    
    %plot(xShift); hold on; plot(zShift); hold off; drawnow;
    
    scatter(xreal,zreal); 
    hold on; 
    scatter(xp,zp); 
    hold off; 
    drawnow;
    
    sse = nanmedian(sqrt((xreal-xp).^2+100*(zreal-zp).^2))
    medist = nanmedian(sqrt((xreal-xp).^2+(zreal-zp).^2))
   % pause(1)
    %sse = nansum(abs(zp-zreal));
end

end