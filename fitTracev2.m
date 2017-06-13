function [est,xshift,zshift] = fitTracev2(xrealArray,zrealArray,truePos,obj)
t = 360-obj.theta;
dmdy = obj.e(2);
xo = obj.opticCentre(1);
zo = obj.opticCentre(2);



st = [0,0,-5,0];
model = @fun;
options = optimset('maxfunevals',2e3*length(st),'maxiter',2e3*length(st),'tolfun',1e-8,'tolx',1e-8);
est = fminsearch(model,st,options);

[~,xshift,zshift] = model(est);


function [sse,xShift,zShift] = fun(params)  
    phi = params(1);
   gamma = params(2);
    %phi = 0;
    dx = params(3);
    dy = params(4);
    
    for k=1:size(xrealArray,1)
        
        xreal = xrealArray(k,:);
        zreal = zrealArray(k,:);

        x = truePos(k,1);
        y = truePos(k,2);
        z = truePos(k,3);
    
        for i=1:length(xreal)
            [xp(i),yp(i),zp(i)] = ARtranslate(x,y,z,-dx,-dy);
            [xp(i),yp(i),zp(i)] = rotArbAxis(xp(i),yp(i),zp(i),phi,gamma,t(i));
            [xp(i),yp(i),zp(i)] = ARtranslate(xp(i),yp(i),zp(i),dx,dy);
            [xp(i),yp(i),zp(i)] = magScale(xp(i),yp(i),zp(i),dmdy.*obj.pxSz/obj.mtz0,xo,zo);
        end

        xreal(or(xp<=-1000,xp>=1000)) = NaN;
        zreal(or(xp<=-1000,zp>=1000)) = NaN;

        xShift(:,k) = xp-xreal;
        zShift(:,k) = zp-zreal;
        
       % scatter3(xreal,yp,zreal,'b'); 
        %hold on; 
        %plot3(xp,yp,zp,'r'); 
        %view([0,0]);
        
    end
    %drawnow;
    %hold off;
    
    %plot(xShift); hold on; plot(zShift); hold off; drawnow;
    
    
%     
    sse = nansum(sqrt((xreal-xp).^2+(zreal-zp).^2));
   % sse = nansum(abs(zShift(:)));
end

end