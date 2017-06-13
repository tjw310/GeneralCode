function recons = test24june(coordsInit,varargin)
%test 24th June
if nargin >3
recons = run(coordsInit,40,-25,varargin);
else
recons = run(coordsInit,40,-25);
end

end


function recons = run(coordsInit,opCentx,opCentz,varargin)
    if nargin>3
        a = varargin{1};
        projections = a{1};
        param = a{2};
        object = a{3};
    else
        param = setParametersv4(opCentx,opCentz);
        [projections,~,~,object]=getBeadConeProjectionsv2(param,coordsInit);
    end
    proj_filtered = filteringv4(projections,param);
    recons = CTbackprojectionv2(proj_filtered, param);
    recons(recons<0)=0;
    %drawmax(recons,param,coordsInit);
    %drawsum(Reconimg,param,coordsInit);
     plotallSum(object,recons,param,coordsInit);
     plotallMax(object,recons,param,coordsInit);    
%     dif = object./max(object(:)) - recons./max(recons(:));
%     drawmax(dif,param,coordsInit);
    drawsliceXY(recons,param,77);
    drawsliceXZ(recons,param,59);
end

function drawsliceXY(recons,param,sliceno)
    figure;
    imagesc(param.rx,param.ry,squeeze(recons(:,:,sliceno)));
    title(sprintf('Slice no = %f',sliceno));
end

function drawsliceXZ(recons,param,sliceno)
    figure;
    imagesc(param.rx,param.ry,squeeze(recons(:,sliceno,:)));
    title(sprintf('Slice no = %f',sliceno));
end

function drawmax(Reconimg,param,coordsInit)
    figure;
    subplot(2,2,1)
    imagesc(param.rx,param.rz,squeeze(max(Reconimg,[],1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); 
    scatter(param.opticCentre(1),param.opticCentre(2),'ro'); hold off;
    subplot(2,2,2)
    imagesc(param.ry,param.rz,squeeze(max(Reconimg,[],2)).'); axis square; xlabel('y'); ylabel('z');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); 
    scatter(0,param.opticCentre(2),'ro'); hold off;
    subplot(2,2,3)
    imagesc(param.rx,param.ry,squeeze(max(Reconimg,[],3))); axis square; xlabel('x'); ylabel('y');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); 
    scatter(param.opticCentre(1),0,'ro'); hold off;
end

function drawsum(Reconimg,param,coordsInit)
    figure;
    subplot(2,2,1)
    imagesc(param.rx,param.rz,squeeze(sum(Reconimg,1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); 
    scatter(param.opticCentre(1),param.opticCentre(2),'ro'); hold off;
    subplot(2,2,2)
    imagesc(param.ry,param.rz,squeeze(sum(Reconimg,2)).'); axis square; xlabel('y'); ylabel('z');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); 
    scatter(0,param.opticCentre(2),'ro'); hold off;
    subplot(2,2,3)
    imagesc(param.rx,param.ry,squeeze(sum(Reconimg,3))); axis square; xlabel('x'); ylabel('y');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); 
    scatter(param.opticCentre(1),0,'ro'); hold off;
end



function plotallSum(object,Reconimg,param,coordsInit)
figure; subplot(2,3,1);
    imagesc(param.rx,param.rz,squeeze(sum(Reconimg,1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); 
    scatter(param.opticCentre(1),param.opticCentre(2),'ro'); hold off;
    subplot(2,3,2)
    imagesc(param.ry,param.rz,squeeze(sum(Reconimg,2)).'); axis square; xlabel('y'); ylabel('z');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); 
    scatter(0,param.opticCentre(2),'ro'); hold off;
    subplot(2,3,3)
    imagesc(param.rx,param.ry,squeeze(sum(Reconimg,3))); axis square; xlabel('x'); ylabel('y');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); 
    scatter(param.opticCentre(1),0,'ro'); hold off;
    subplot(2,3,4);
    imagesc(param.rx,param.rz,squeeze(sum(object,1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); 
    scatter(param.opticCentre(1),param.opticCentre(2),'ro'); hold off;
    subplot(2,3,5)
    imagesc(param.ry,param.rz,squeeze(sum(object,2)).'); axis square; xlabel('y'); ylabel('z');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); 
    scatter(0,param.opticCentre(2),'ro'); hold off;
    subplot(2,3,6)
    imagesc(param.rx,param.ry,squeeze(sum(object,3))); axis square; xlabel('x'); ylabel('y');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); 
    scatter(param.opticCentre(1),0,'ro'); hold off;
end

function plotallMax(object,Reconimg,param,coordsInit)
figure; subplot(2,3,1);
    imagesc(param.rx,param.rz,squeeze(max(Reconimg,[],1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); 
    scatter(param.opticCentre(1),param.opticCentre(2),'ro'); hold off;
    subplot(2,3,2)
    imagesc(param.ry,param.rz,squeeze(max(Reconimg,[],2)).'); axis square; xlabel('y'); ylabel('z');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); 
    scatter(0,param.opticCentre(2),'ro'); hold off;
    subplot(2,3,3)
    imagesc(param.rx,param.ry,squeeze(max(Reconimg,[],3))); axis square; xlabel('x'); ylabel('y');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); 
    scatter(param.opticCentre(1),0,'ro'); hold off;
    subplot(2,3,4);
    imagesc(param.rx,param.rz,squeeze(max(object,[],1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); 
    scatter(param.opticCentre(1),param.opticCentre(2),'ro'); hold off;
    subplot(2,3,5)
    imagesc(param.ry,param.rz,squeeze(max(object,[],2)).'); axis square; xlabel('y'); ylabel('z');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); 
    scatter(0,param.opticCentre(2),'ro'); hold off;
    subplot(2,3,6)
    imagesc(param.rx,param.ry,squeeze(max(object,[],3))); axis square; xlabel('x'); ylabel('y');hold on;
    set(gca,'ydir','normal');
    scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); 
    scatter(param.opticCentre(1),0,'ro'); hold off;
end
