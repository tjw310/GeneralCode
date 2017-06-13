function [projections,tracks,coordsInit,object] = getBeadConeProjectionsv2(param,varargin)
if nargin>1
    bool = 1;
    coordsInit = varargin{1};
else
    bool=0;
end
projections = [];

%x,y,z coords in object space, real size.
sc = single((-param.nPx/2:param.nPx/2-1)*param.pxSz/param.mtz0);
dpx = param.pxSz/param.mtz0;
[x,z] = meshgrid(single(sc));
z = gpuArray(z);
x = gpuArray(x);
y = gpuArray(single(sc));

%Bead coords, IN this case the object rotates about the z axis. Depth is y
nBeads = 10;
rBeads = 3.5e-3; %real radius of beads in mm

nslices = ceil(2*rBeads/dpx);
    
tracks =[];
figure;
for k=1:length(param.theta);
    disp(k/length(param.theta)*100);
    if k==1
        if bool~=1
        coordsInit = getInitialObject(nBeads,param.w); %get object initial coords
        end
        coords = coordsInit;
    else
        coords = rotateCoords(coordsInit,param.theta(k));
    end
    
    %create object space and fill it will beads at coords
    
    
    allbeads = gpuArray(single(zeros(param.nPx,param.nPx)));
    XY = gpuArray(single(zeros(param.nPx,param.nPx)));
    YZ = gpuArray(single(zeros(param.nPx,param.nPx)));
    for i=1:nBeads
        yidx = coords(2,i);
        [~,lc]=min(abs(y-yidx));
        out = gpuArray(single(zeros(param.nPx,param.nPx)));
        XYbeads = gpuArray(single(zeros(param.nPx,param.nPx)));
        YZbeads = gpuArray(single(zeros(param.nPx,param.nPx)));
        N=0;
        for j=-ceil(nslices/2):1:ceil(nslices/2)
            object = gpuArray(single(zeros(param.nPx,param.nPx)));
            n1=0;
            n2=0;
            if and(j+lc>0,j+lc<length(sc))
                yC = y(j+lc);
                [xidx,yC,zidx]=coordsMagChange(coords(1,i),yC,coords(3,i),param);     
            r2 = sqrt(rBeads^2-(dpx*abs(j))^2)*(1+param.dmdy*yC);
            r = sqrt((xidx-x).^2+...
                    (zidx-z).^2);
                if nnz(r<=r2)>0
                    object(r<=r2)=1;
                    n1 = nnz(r<=r2);
                    if nnz(and(r>r2,(r-r2)<dpx))>0
                        object(and(r>r2,(r-r2)<dpx)) = (1-(r(and(r>r2,(r-r2)<dpx))-r2)./dpx);
                        n2 = sum(object(:));
                    end
                    out = out+object;
                end
            end
            N = N +n1+n2;
            XYbeads(j+lc,:) = sum(object,1);
            YZbeads(:,j+lc) = sum(object,2);
        end
         XY = XY + XYbeads./N;
         YZ = YZ + YZbeads./N;
         allbeads = allbeads + out./N;    
    end

    %records tracks
    tracks(1,:,k) = scalex(coords(1,:),coords(2,:),param);
    tracks(2,:,k) = coords(2,:);
    tracks(3,:,k) = scalez(coords(3,:),coords(2,:),param);
    
    %get xz projection
    projections(:,:,k) = gather(allbeads.');
    
    
    
    subplot(2,2,1)
    imagesc(sc,sc,allbeads); xlabel('x /mm');ylabel('z /mm'); axis square;
    set(gca,'ydir','normal');
    hold on
    scatter(param.opticCentre(1),param.opticCentre(2),'ro');
    scatter(scalex(coords(1,:),coords(2,:),param),scalez(coords(3,:),coords(2,:),param),'rx');
    hold off
    subplot(2,2,2)
    imagesc(sc,sc,XY); axis square; set(gca,'ydir','normal'); xlabel('x /mm');ylabel('y /mm');
    hold on
    scatter(param.opticCentre(1),0,'ro');
    scatter(scalex(coords(1,:),coords(2,:),param),coords(2,:),'rx');
    hold off
    subplot(2,2,3)
    imagesc(sc,sc,YZ); axis square; set(gca,'ydir','normal'); xlabel('y /mm');ylabel('z /mm');
    hold on
    scatter(0,param.opticCentre(2),'ro');
    scatter(coords(2,:),scalez(coords(3,:),coords(2,:),param),'rx');
    hold off
    drawnow;
end

[x2,y2,z2]=meshgrid(single(sc));
    object = gpuArray(single(zeros(param.nPx,param.nPx,param.nPx)));
for i=1:nBeads
        xC = coordsInit(1,i);
        zC = coordsInit(3,i);
        yC = coordsInit(2,i);
        r = sqrt((xC-x2).^2+...
                (yC-y2).^2+(zC-z2).^2);
            if nnz(r<=rBeads)>0
                object(r<=rBeads)=1;
                if nnz(and(r>rBeads,(r-rBeads)<dpx))>0
                    object(and(r>rBeads,(r-rBeads)<dpx)) = (1-(r(and(r>rBeads,(r-rBeads)<dpx))-rBeads)./dpx);
                end
            end
end
    figure;
    imagesc(sc,sc,(squeeze(sum(object,1))).'); axis square; xlabel('x /mm'); ylabel('y /mm');
    set(gca,'ydir','normal'); hold on;
    scatter(param.opticCentre(1),param.opticCentre(2),'ro');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx');
    hold off
    set(gca,'ydir','normal');
    
    figure;
    for l=1:size(tracks,2)   
        subplot(1,3,3)
        plot(squeeze(tracks(1,l,:))/param.w*param.nPx,squeeze(tracks(2,l,:))/param.w*param.nPx); hold on;
        if l==1
            scatter(param.opticCentreP(1),0,'rx');
        end
        xlabel('x');ylabel('y');
        subplot(1,3,2)
        plot(squeeze(tracks(2,l,:))/param.w*param.nPx,squeeze(tracks(3,l,:))/param.w*param.nPx); hold on;
        if l==1
            scatter(0,param.opticCentreP(2),'rx');
        end
        xlabel('y'); ylabel('z');
        subplot(1,3,1)
        plot(squeeze(tracks(1,l,:))/param.w*param.nPx,squeeze(tracks(3,l,:))/param.w*param.nPx); hold on;
        if l==1
            scatter(param.opticCentreP(1),param.opticCentreP(2),'rx');
        end
        xlabel('x'); ylabel('z');
    end
    hold off;
    
    figure;
    scatter3(coordsInit(1,:),coordsInit(2,:),coordsInit(3,:),'filled');
    xlim([min([param.rx,sc]),max([param.rx,sc])]);
    ylim([min([param.ry,sc]),max([param.ry,sc])]);
    zlim([min([param.rz,sc]),max([param.rz,sc])]);
    xlabel('x'); ylabel('y'); zlabel('z');
end

function coords = getInitialObject(nBeads,w)
    for i=1:nBeads
        %creates beads randomly (i.e can overlap);
        coords(1,i) = (2*rand(1)-1)*w/2/sqrt(2)*0.9;
        coords(2,i) = (2*rand(1)-1)*w/2/sqrt(2)*0.9;
        coords(3,i) = (2*rand(1)-1)*w/2/sqrt(2)*0.9;
    end
end

function out = scalex(x,y,param)
    out = (x-param.opticCentre(1)).*(1+param.dmdy*y)+param.opticCentre(1);
end

function out = scalez(z,y,param)
    out = (z-param.opticCentre(2)).*(1+param.dmdy*y)+param.opticCentre(2);
end

function [x,y,z] = coordsMagChange(x,y,z,param)
        %find depth position and calculates relative change in magnification
        %and then scales the x and z positions respectively.
        dm = 1+param.dmdy*y;
        x = (x-param.opticCentre(1)).*dm+param.opticCentre(1);
        z = (z-param.opticCentre(2)).*dm+param.opticCentre(2);
end

function testMagChange(coords)
    %test magnifciation change
    figure;
    for j=1:5
        dmrdy = 0.2+j*0.2;
        dm = 1+coords(2,:)*dmrdy;
        coords(1,:) = (coords(1,:)-param.opticCentre(1)).*dm+param.opticCentre(1);
        coords(3,:) = (coords(3,:)-param.opticCentre(2)).*dm+param.opticCentre(2);
        scatter3(coords(1,:),coords(2,:),coords(3,:),'filled');
        hold on
        xlabel('x');
        ylabel('y');
        zlabel('z');
        drawnow;
    end
    scatter3(param.opticCentre(1),0,param.opticCentre(2),'rx');
    hold off
    axis square
end

function coordsOut = rotateCoords(coords,theta)
%rotates coords about the z-axis (so assumes AoR is z-axis). Need to
%convert beta to radians
theta = theta/180*pi;
coordsOut(1,:) = coords(1,:).*cos(theta)-coords(2,:).*sin(theta);
coordsOut(2,:) = coords(1,:).*sin(theta)+coords(2,:).*cos(theta);
coordsOut(3,:) = coords(3,:);
end
