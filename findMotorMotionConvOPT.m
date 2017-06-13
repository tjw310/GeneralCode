function [sse,traces,motorMotion,AoRoffset,zMotion,zPos,AoRangle] = findMotorMotionConvOPT(obj)

[traces,zMot] = getTraces(obj);

start_guesses = [obj.mtz0,0,0];
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
theta = linspace(0,2*pi-2*pi/obj.nProj,obj.nProj);

count = 1;
for k=1:size(traces,1)
    xdata = traces(k,:);
    estimates = fminsearch(model,start_guesses,options);
    [sse(k),motMot,fit] = model(estimates);
    
    if sse(k)<1000 && sse(k)>1
        fits(count,:) = fit;
        xdat(count,:) = xdata;
        AoRoffsets(count) = estimates(3);
        motorMotion(count,:) = motMot;
        zMotion(count,:) = zMot(k,:)-nanmean(zMot(k,:));
        zPos(count) = nanmean(zMot(k,:));
        count = count+1;
    end
end

mip = squeeze(max(obj.projections,[],3));
fc = 0.2;
mip = (mip-min(mip(:)))./(max(mip(:))-min(mip(:)))/fc;
mip(mip>1) = 1;
mip = mip.*(max(zPos)-min(zPos))+min(zPos);
figure; subplot(1,2,1); imagesc(obj.x,obj.z,mip.'); set(gca,'ydir','normal'); axis square; 
subplot(1,2,2); imagesc(obj.x,obj.z,mip.'); set(gca,'ydir','normal');  axis square; hold on;
for l=1:size(traces,1)
   if sse(l)<1000 && sse(l)>1
    scatter(traces(l,:),zMot(l,:),[],repmat(zPos(l),1,length(traces(l,:))),'filled','markeredgecolor','r'); 
   end
end
scatter(xdat(:,1),zPos.'+zMotion(:,1),[],zPos,'d','filled','markeredgecolor','k');
hold off;
savefig(gcf,strcat(obj.path,sprintf('mipandtraces.fig')));

e = fitLinear(zPos.',AoRoffsets.'); 
AoRangle = atan(e(1))/pi*180;
AoRoffset = e(2);
xlabel('AoR X-Offset Position /pixel'); ylabel('Microsphere Z Position /pixel');
title(sprintf('AoRangle=%0.3f ,AoRoffset=%.1f',AoRangle,AoRoffset));
savefig(gcf,strcat(obj.path,sprintf('/AoRangle=%0.3f ,AoRoffset=%.1f.fig',AoRangle,AoRoffset)));

figure;
zPosNorm = (zPos-min(zPos))./(max(zPos)-min(zPos));
cmap = colormap;
for p=1:size(motorMotion,1)
    [~,lc] = min(abs(zPosNorm(p)-linspace(0,1,size(cmap,1))));
    plot(motorMotion(p,:),linspace(0,360,obj.nProj),'color',cmap(lc,:)); hold on;
end
hold off; xlabel('Microsphere Excess X-Motion /pixel'); ylabel('Projection Angle');
savefig(gcf,strcat(obj.path,'/motorMotion.fig'));

figure;
for p=1:size(zMotion,1)
    [~,lc] = min(abs(zPosNorm(p)-linspace(0,1,size(cmap,1))));
    plot(zMotion(p,:),linspace(0,360,obj.nProj),'color',cmap(lc,:)); hold on;
end
hold off; xlabel('Microsphere Z-Motion /pixel'); ylabel('Projection Angle');
savefig(gcf,strcat(obj.path,'/zMotion.fig'));

figure;
for p=1:size(fits,1)
    [~,lc] = min(abs(zPosNorm(p)-linspace(0,1,size(cmap,1))));
    scatter(xdat(p,:),linspace(0,360,obj.nProj),[],cmap(lc,:)); hold on;
    plot(fits(p,:),linspace(0,360,obj.nProj),'r');
end
hold off; xlabel('Microsphere X-Motion /pixel'); ylabel('Projection Angle');
savefig(gcf,strcat(obj.path,'/xMotion.fig'));


    function [sse,motMot,fit] = fun(params)
            A = params(1);
            B = params(2);
            C = params(3);
            fit = A.*sin(theta+B)+C;
            %scatter(1:400,xdata); hold on; plot(fit); hold off; drawnow;
            sse = nansum((xdata-fit).^2);
            motMot = xdata-fit;
    end


    function [x,z,bool] = getTraces(obj)
        if isempty(obj.peaksX)
            [x,z,~,~] = obj.findBeadsAllSlices;
        else
        x = obj.peaksX; z = obj.peaksZ;
        end

        bool(1:size(x,1),1) = ones(1,size(x,1));
        for i=1:size(x,2)-1
            [xpair,zpair] = findClosestPairs(x(:,i),z(:,i),x(:,i+1),z(:,i+1));

            x(:,i+1) = xpair; z(:,i+1) = zpair;
            
            bool(:,i+1) = ~isnan(x(:,i+1));            
            x(isnan(x(:,i+1)),i+1) = x(isnan(x(:,i+1)),i);
            z(isnan(z(:,i+1)),i+1) = z(isnan(z(:,i+1)),i);
        end
        
        x(bool==0) = NaN;
        z(bool==0) = NaN;
    end
            
end