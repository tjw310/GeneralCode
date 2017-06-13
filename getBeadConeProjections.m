function [projections,tracks,coordsInit] = getBeadConeProjections(param,varargin)
if nargin>1
    bool = 1;
    coordsInit = varargin{1};
else
    bool=0;
end
projections = [];

%x,y,z coords in object space, real size.
sc = single(linspace(-param.w/2,param.w/2,param.nPx));
dpx = abs(sc(2)-sc(1));
[x,y,z] = meshgrid(sc);
y = gpuArray(y);
x = gpuArray(x);
z = gpuArray(z);

%Bead coords, IN this case the object rotates about the z axis. Depth is y
nBeads = 10;
rBeads = 1e-3; %real radius of beads in mm

tracks =[];
figure;
for k=1:length(param.theta);
    if k==1
        if bool~=1
        coordsInit = getInitialObject(nBeads,param.w); %get object initial coords
        end
        coords = coordsMagChange(coordsInit,param); %scale coords due to magnification
    else
        coords = coordsMagChange(rotateCoords(coordsInit,param.theta(k)),param);
    end
    
    %create object space and fill it will beads at coords
    ob = gpuArray(zeros(param.nPx,param.nPx,param.nPx));
    for i=1:nBeads
        r = sqrt((((coords(1,i)-x+param.opticCentre(1)).*(1+param.dmdy*coords(2,i)))+param.opticCentre(1)).^2+...
            (coords(2,i)-y).^2+...
                (((coords(3,i)-z+param.opticCentre(2)).*(1+param.dmdy*coords(2,i)))+param.opticCentre(2)).^2);
        r = sqrt(((coords(1,i)-x).^2+...
        (coords(2,i)-y).^2+...
            ((coords(3,i)-z)).^2));
        n1 = nnz(r<=rBeads);
        dif = (r-rBeads)./dpx;
        n2 = sum(dif(and(dif>0,dif<1)));
        tot = n1+n2;
        ob(dif<=0)=1/tot;
        %ob(and(dif>0,dif<1)) = dif(and(dif>0,dif<1))./tot;
        
%         [~,lc] = min(abs(sc-coords(2,i)));
%         imagesc(squeeze(ob(lc,:,:)));
%         drawnow;
%         pause(2);
%         
   
    end

    %records tracks
    tracks(:,:,k) = coords;
    
    %get xz projection
    projections(:,:,k) = gather(squeeze(sum(ob,1)));

    imagesc(sc,sc,squeeze(sum(ob,1)).'); xlabel('x /mm');ylabel('y /mm');
    set(gca,'ydir','normal');
    hold on
    scatter(param.opticCentre(1),param.opticCentre(2),'ro');
    scatter(coords(1,:),coords(3,:),'rx');
    hold off
    drawnow;
end
    
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
    xlabel('x'); ylabel('y'); zlabel('z');

end

function coords = getInitialObject(nBeads,w)
    for i=1:nBeads
        %creates beads randomly (i.e can overlap);
        coords(1,i) = (2*rand(1)-1)*w/2/sqrt(2);
        coords(2,i) = (2*rand(1)-1)*w/2/sqrt(2);
        coords(3,i) = (2*rand(1)-1)*w/2/sqrt(2);
    end
end

function coordsOut = coordsMagChange(coords,param)
     for i=1:size(coords,2)
        %find depth position and calculates relative change in magnification
        %and then scales the x and z positions respectively.
        dm = 1+param.dmdy*coords(2,i);
        coordsOut(1,i) = (coords(1,i)-param.opticCentre(1)).*dm+param.opticCentre(1);
        coordsOut(3,i) = (coords(3,i)-param.opticCentre(2)).*dm+param.opticCentre(2);
        coordsOut(2,i) = coords(2,i);
     end
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
