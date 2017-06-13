function [ img ] = CTbackprojectionv3( proj, param )
%CTBACKPROJECTION adapated from mathworks. Optimised for GPU.

img = zeros(param.nPx, param.nPx, param.nPx, 'single');

[x,y] = meshgrid((-param.nPx/2:param.nPx/2-1+param.opticCentreP(1)),(-param.nPx/2:param.nPx/2-1));
if param.gpu == 1
        x = gpuArray(x);
    y = gpuArray(y);
end
figure;
for i = 1:param.nProj
    disp(i);
    tic
   img = img + backprojection(single(proj(:,:,i)),param,i,x,y);
   t(i) = toc;
end

figure;
plot(t);
title(strcat(num2str(mean(t)),'\pm',num2str(std(t))));

end

function vol = backprojection(proj,param,iview,x,y)

beta = param.theta(iview)/360*2*pi;
vol = single(zeros(param.nPx,param.nPx,param.nPx,'single'));

u = single(param.D.*(x*cos(beta)+y*sin(beta)+sqrt(param.R^2-param.D^2))...
    ./(x*sin(beta)-y*cos(beta)+param.D));

u = u(:,1:param.nPx);
frac = single(param.D./(x*sin(beta)-y*cos(beta)+param.D));
frac = frac(:,1:param.nPx);

if param.gpu == 1
    proj = gpuArray(proj);
    ys = gpuArray(param.y(1)-1);
    u2 = gpuArray(u-param.x(1)+1);
end

for iz = 1:param.nPx   
    if param.gpu == 1
        v = frac.*(param.z(iz))+param.opticCentreP(2);
        v2 = v-ys;
        %vol(:,:,iz) = gather(param.D^2.*frac./(param.D^2+v.^2).*interp2(proj',circshift(u2,[0,param.opticCentreP(1)]),v2,param.interptype));
        vol(:,:,iz) = flipud(gather(param.D^2.*frac./(param.D^2+v.^2).*interp2(proj',u2,v2,param.interptype)));
    else
        v = single(frac.*param.z(iz)+param.off_v);
        vol(:,:,iz) = param.D^2.*frac./(param.D^2+v.^2).*interp2(proj',u-param.x(1)+1,v-param.y(1)+1,param.interptype);
    end 

    
end
% subplot(1,4,1)
% imagesc(u2);
% colorbar;
% subplot(1,4,2);
% imagesc(circshift(u2,[0,param.opticCentreP(1)]));
% colorbar;
% subplot(1,4,3);
% imagesc(abs(u2-circshift(u2,[0,param.opticCentreP(1)])));
% drawnow;
% u3 = circshift(u2,[0,param.opticCentreP(1)]);
% subplot(1,4,4);
% plot(u3(:,50)-u3(:,51))

vol(isnan(vol))=0;

return

end

