function estimates = fit180projectionsv2(obj,fixed,X,Z,traceEst,angleEst)
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
    if mn<50
        yinit(k) = yf(lc);
        xinit(k) = X(i,1);
        zinit(k) = Z(i,1);
        k=k+1;
    end
    end
end
figure; scatter(xf,zf); hold on; scatter(xinit,zinit,'rx'); hold off; drawnow;
figure;
magt = traceEst(2);
st = [traceEst(1),traceEst(3),traceEst(4),angleEst(1),angleEst(1),magt];
model = @fun;

options = optimset('maxfuneval',2e4,'maxiter',2e4,'tolx',1e-6,'tolfun',1e-6);

%[sse,medDist,plotx,plotz] = fun(st);

estimates = fminsearch(model,st,options)
[sse,medDist,plotx,plotz] = model(estimates);
figure;

for t=1:size(X,1)
    scatter(X(t,:),Z(t,:),[],1:400,'filled'); hold on;
    scatter(X(t,:),Z(t,:),'b'); hold on;
end

%scatter(plotx,plotz,[],1:400,'filled');
scatter(plotx,plotz,'rx');
%ylim([-900,-750]);xlim([-850,850]);
title(num2str(medDist));
hold off; drawnow;

function [sse,medDist,plotx,plotz] = fun(params) %common global parameters for the fit.    
    xo = obj.opticCentre(1);
    zo = obj.opticCentre(2);
    
    dmdy = params(1)*6.5e-3/params(6);
    alpha = params(4);
    beta = params(5);
    
    xa = params(2);
    ya = params(3);
    mt = params(6);

    
    [u,v,w] = getAxDrVc(alpha,beta);

    Xcopy = X;
    Zcopy = Z;

    model2 = @fitTrace;

     plotx = []; plotz = [];
     for t=1:size(X,1)
    scatter(X(t,:),Z(t,:),[],1:400,'filled'); hold on;
    scatter(X(t,:),Z(t,:),'b'); hold on;
    end

    for beadNo = 1:length(yinit)
%        % figure;
        st2 = [yinit(beadNo)/6.5e-3*mt,xinit(beadNo),zinit(beadNo)];
         estimates2 = fminsearch(model2,st2);
         [tracesse(beadNo),x2pred,z2pred,medDist] = model2(estimates2);
         
         %[tracesse(beadNo),x2pred,z2pred,medDist] = fitTrace(st2);
         plotx = horzcat(plotx,x2pred);
         plotz = horzcat(plotz,z2pred);

%        hold on; 
%        scatter(x2pred,z2pred,[],1:400,'filled');
%        scatter(x2pred,z2pred,'rx');
%        scatter(x2pred(1),z2pred(1),'g','filled'); 
%        scatter(xo,zo,'k','filled'); hold off
%         
%        %ylim([-625,-580]);xlim([-1050,1050]);
%        
%        drawnow;
      % pause(10);
    end
    
    sprintf('dmdy=%.3f,xa=%.3f,ya=%.3f,alpha=%.3f,beta=%.3f,mt=%.3f',params(1),params(2),params(3),params(4),params(5),params(6))
    

    
   %tracesse =0;
    sse = nanmean(tracesse)

    function [sse2,x2pred,z2pred,medDist] = fitTrace(params) %parameters specific to each bead (ylocation)
        theta = (360-obj.theta)/180*pi;

        y = params(1);
        
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
        
       % plot(theta); 
       
%        if beadNo==1
%        for t=1:67
%             scatter(Xcopy(t,:),Zcopy(t,:),[],1:400,'filled'); hold on;
%             scatter(Xcopy(t,:),Zcopy(t,:),'b'); hold on;
%        end
%        hold on; 
%        scatter(x2pred,z2pred,[],1:400,'filled');
%        scatter(x2pred,z2pred,'rx');
%        scatter(x2pred(1),z2pred(1),'g','filled'); 
%        scatter(xo,zo,'k','filled'); hold off
%        ylim([-900,-750]);xlim([-850,850]);
%        drawnow;
%        end
        %imagesc(r2array); drawnow; 
       %pause(10);

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
        
        %scatter(xclose,zclose); hold on; scatter(x2pred,z2pred,'rx'); hold off;
        %drawnow;
        sse2 = nanmean(mn2);
        %sse2 = nanmedian(mn2);
        medDist = nanmedian(mn2);
    end
end
end
