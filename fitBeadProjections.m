function [traceEstimates,angleEst,dx,dz,index,yinit] = fitBeadProjections(obj,fixed,X,Z,traceGuess,angleGuess)
%[X,Z] = findBeadsAllSlices(obj);

X(or(X<-1018,X>1018)) = NaN;
Z(or(Z<-1018,Z>1018)) = NaN;

for j=1:length(fixed)
    xf(j) = fixed(j).centreX-2560/2;
    zf(j) = fixed(j).centreY-2160/2;
    yf(j) = fixed(j).zDepth;
end


k=1;
for i=1:length(X(:,1))
    if ~isnan(X(i,1)) && ~isnan(Z(i,1))
     r = sqrt((xf-X(i,1)).^2+(zf-Z(i,1)).^2);
    [mn,lc] = min(r);
    if mn<150
        yinit(k) = yf(lc);
        xinit(k) = X(i,1);
        zinit(k) = Z(i,1);
        index(k) = i;
        k=k+1;
    end
    end
end

traceEstimates = traceGuess;
angleEst = angleGuess;
modelMotion = @fitMotion;
modelXYZ = @fitXYZ;
modelAngle= @fitAngles;
modelY = @fitY;


options = optimset('maxfuneval',2e4,'maxiter',2e4,'tolx',1e-6,'tolfun',1e-6);
figure;
for l=1
    
%for i=1:length(xinit)
for i=6;
    beadNo = i;
    %if l==1
    %fitY(yinit(i));
       posEst = fminsearch(modelY,yinit(i),options);
    %else
      %  posEst = fminsearch(modelXYZ,[yinit(i),xinit(i),zinit(i)],options);     
      %  xinit(i) = posEst(2);
      %  zinit(i) = posEst(3);
    %end
    yinit(i) = posEst(1);
end

traceEstimates = fminsearch(modelMotion,traceEstimates,options)
angleEst = fminsearch(modelAngle,angleEst,options)
[~,dx,dz,xplot,zplot] = modelAngle(angleEst);
end

figure;
for t=1:size(X,1)
    scatter(X(t,:),Z(t,:),[],1:400,'filled'); hold on;
    scatter(X(t,:),Z(t,:),'b'); hold on;
end

scatter(xplot,zplot,'rx');
%ylim([-900,-750]);xlim([-850,850]);
%title(num2str(medDist));
hold off; drawnow;

function [sse,medDist] = fitMotion(params) %common global parameters for the fit.    
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = params(1)*6.5e-3/params(2);
    alpha = angleEst(1);
    beta = angleEst(2);
    
    xa = params(3);
    ya = params(4);
    mt = params(2);
  
    [u,v,w] = getAxDrVc(alpha,beta);

    Xcopy = X;
    Zcopy = Z;
    
    %theta = (360-obj.theta)/180*pi;
    theta = obj.theta/180*pi;
    for beadNo = 1:length(xinit)
        y = yinit(beadNo)/6.5e-3*mt;
        x = (xinit(beadNo)-xo).*(1-dmdy*y)+xo;
        z = (zinit(beadNo)-zo).*(1-dmdy*y)+zo;

        f = (u*(x-xa)+v*(y-ya)+z).*(1-cos(theta));

        rx = u*f+(x-xa).*cos(theta)+(-w*(y-ya)+v*z).*sin(theta)+xa;
        ry = v*f+(y-ya).*cos(theta)+(w*(x-xa)-u*z).*sin(theta)+ya;
        rz = w*f+z.*cos(theta)+(-v*(x-xa)+u*(y-ya)).*sin(theta);

        x2pred = (rx-xo)./(1-dmdy*ry)+xo;
        z2pred = (rz-zo)./(1-dmdy*ry)+zo;

        x2array = repmat(x2pred,size(Xcopy,1),1);
        z2array = repmat(z2pred,size(Zcopy,1),1);
              
        r2array = sqrt((x2array-Xcopy).^2+(z2array-Zcopy).^2);
        
        for n=1:size(r2array,2)
            [mn2(n),lc2] = min(r2array(:,n));
            xclose(n) = Xcopy(lc2,n);
            zclose(n) = Zcopy(lc2,n);
              %line = Xcopy(:,n);
            %lfill = line((1:length(line))~=lc2);
            %Xcopy(:,n) = lfill;
            %line = Zcopy(:,n);
            %lfill = line((1:length(line))~=lc2);
            %Zcopy(:,n) = lfill;    
        end
        
       if beadNo==6
           scatter(xclose,zclose,'bo'); hold on
           scatter(xclose,zclose,[],1:400,'filled');
           scatter(x2pred,z2pred,[],1:400,'filled');
           scatter(x2pred,z2pred,'rx'); 
           scatter(x2pred(1),z2pred(1),'g','filled');
           hold off;
           xlim([min(x2pred)-100,max(x2pred)+100]);
           ylim([min(z2pred)-100,max(z2pred)+100]);
           title('Find Motion');
           drawnow;
       end

        err(beadNo) = nanmedian(mn2);
    end

    sse = nanmedian(err);
      
end

function [sse,medDist] = fitY(params) %common global parameters for the fit.    
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = traceEstimates(1)*6.5e-3/traceEstimates(2);
    alpha = angleEst(1);
    beta = angleEst(2);
    
    xa = traceEstimates(3);
    ya = traceEstimates(4);
    mt = traceEstimates(2);
  
    [u,v,w] = getAxDrVc(alpha,beta);

    Xcopy = X;
    Zcopy = Z;
    
    %theta = (360-obj.theta)/180*pi;
    theta = obj.theta/180*pi;
    y = params(1)/6.5e-3*mt;
    x = (xinit(beadNo)-xo).*(1-dmdy*y)+xo;
    z = (zinit(beadNo)-zo).*(1-dmdy*y)+zo;

    f = (u*(x-xa)+v*(y-ya)+z).*(1-cos(theta));

    rx = u*f+(x-xa).*cos(theta)+(-w*(y-ya)+v*z).*sin(theta)+xa;
    ry = v*f+(y-ya).*cos(theta)+(w*(x-xa)-u*z).*sin(theta)+ya;
    rz = w*f+z.*cos(theta)+(-v*(x-xa)+u*(y-ya)).*sin(theta);

    x2pred = (rx-xo)./(1-dmdy*ry)+xo;
    z2pred = (rz-zo)./(1-dmdy*ry)+zo;
    


    x2array = repmat(x2pred,size(Xcopy,1),1);
    z2array = repmat(z2pred,size(Zcopy,1),1);

    r2array = sqrt((x2array-Xcopy).^2+25.*(z2array-Zcopy).^2); 

    for n=1:size(r2array,2)
        [mn2(n),lc2] = min(r2array(:,n));
        xclose(n) = Xcopy(lc2,n);
        zclose(n) = Zcopy(lc2,n);  
    end
    
        %if beadNo==7
           scatter(xclose,zclose,'bo'); hold on
           scatter(xclose,zclose,[],1:400,'filled');
           scatter(x2pred,z2pred,[],1:400,'filled');
           scatter(x2pred,z2pred,'rx'); 
           scatter(x2pred(1),z2pred(1),'g','filled');
           hold off;
           xlim([min(x2pred)-100,max(x2pred)+100]);
           ylim([min(z2pred)-100,max(z2pred)+100]);
           title('Find Best Y Pos');
           drawnow;
      % end

    sse = nanmedian(mn2);
end

function [sse,dx,dz,plotx,plotz] = fitAngles(params) %common global parameters for the fit.    
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = traceEstimates(1)*6.5e-3/traceEstimates(2);
    alpha = params(1);
    beta = params(2);
    
    xa = traceEstimates(3);
    ya = traceEstimates(4);
    mt = traceEstimates(2);
  
    [u,v,w] = getAxDrVc(alpha,beta);

    Xcopy = X;
    Zcopy = Z;
    plotx = []; plotz=[];
    %theta = (360-obj.theta)/180*pi;
    theta = obj.theta/180*pi;
    for beadNo = 1:length(xinit)
    %for beadNo = 5
        y = yinit(beadNo)/6.5e-3*mt;
        x = (xinit(beadNo)-xo).*(1-dmdy*y)+xo;
        z = (zinit(beadNo)-zo).*(1-dmdy*y)+zo;

        f = (u*(x-xa)+v*(y-ya)+z).*(1-cos(theta));

        rx = u*f+(x-xa).*cos(theta)+(-w*(y-ya)+v*z).*sin(theta)+xa;
        ry = v*f+(y-ya).*cos(theta)+(w*(x-xa)-u*z).*sin(theta)+ya;
        rz = w*f+z.*cos(theta)+(-v*(x-xa)+u*(y-ya)).*sin(theta);

        x2pred = (rx-xo)./(1-dmdy*ry)+xo;
        z2pred = (rz-zo)./(1-dmdy*ry)+zo;

        x2array = repmat(x2pred,size(Xcopy,1),1);
        z2array = repmat(z2pred,size(Zcopy,1),1);
        
        plotx = horzcat(plotx,x2pred);
        plotz = horzcat(plotz,z2pred);      
       
        r2array = sqrt((x2array-Xcopy).^2+(z2array-Zcopy).^2);
        
        for n=1:size(r2array,2)
            [mn2(n),lc2] = min(r2array(:,n));
            xclose(n) = Xcopy(lc2,n);
            zclose(n) = Zcopy(lc2,n);
            dx(beadNo,n) = Xcopy(lc2,n)-x2array(lc2,n);
            dz(beadNo,n) = Zcopy(lc2,n)-z2array(lc2,n);
              %line = Xcopy(:,n);
            %lfill = line((1:length(line))~=lc2);
            %Xcopy(:,n) = lfill;
            %line = Zcopy(:,n);
            %lfill = line((1:length(line))~=lc2);
            %Zcopy(:,n) = lfill;    
        end
        
       if beadNo==6
           scatter(xclose,zclose,'bo'); hold on
           scatter(xclose,zclose,[],1:400,'filled');
           scatter(x2pred,z2pred,[],1:400,'filled');
           scatter(x2pred,z2pred,'rx'); 
           scatter(x2pred(1),z2pred(1),'g','filled');
           hold off;
           xlim([min(x2pred)-100,max(x2pred)+100]);
           ylim([min(z2pred)-100,max(z2pred)+100]);
           title('Find Angles');
           drawnow;
       end
 
        err(beadNo) = nanmedian(mn2);
    end
    
    sse = nanmedian(err);
      
end

function [sse,medDist] = fitXYZ(params) %common global parameters for the fit.    
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = traceEstimates(1)*6.5e-3/traceEstimates(2);
    alpha = angleEst(1);
    beta = angleEst(2);
    
    xa = traceEstimates(3);
    ya = traceEstimates(4);
    mt = traceEstimates(2);
  
    [u,v,w] = getAxDrVc(alpha,beta);

    Xcopy = X;
    Zcopy = Z;
    
    theta = (360-obj.theta)/180*pi;
    y = params(1)/6.5e-3*mt;
    x = (params(2)-xo).*(1-dmdy*y)+xo;
    z = (params(3)-zo).*(1-dmdy*y)+zo;

    f = (u*(x-xa)+v*(y-ya)+z).*(1-cos(theta));

    rx = u*f+(x-xa).*cos(theta)+(-w*(y-ya)+v*z).*sin(theta)+xa;
    ry = v*f+(y-ya).*cos(theta)+(w*(x-xa)-u*z).*sin(theta)+ya;
    rz = w*f+z.*cos(theta)+(-v*(x-xa)+u*(y-ya)).*sin(theta);

    x2pred = (rx-xo)./(1-dmdy*ry)+xo;
    z2pred = (rz-zo)./(1-dmdy*ry)+zo;

    x2array = repmat(x2pred,size(Xcopy,1),1);
    z2array = repmat(z2pred,size(Zcopy,1),1);

    r2array = sqrt((x2array-Xcopy).^2+(z2array-Zcopy).^2);

    for n=1:size(r2array,2)
        [mn2(n),lc2] = min(r2array(:,n));
        xclose(n) = Xcopy(lc2,n);
        zclose(n) = Zcopy(lc2,n);  
    end

    sse = nanmedian(mn2);
end
end