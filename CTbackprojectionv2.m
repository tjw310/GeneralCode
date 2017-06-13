function [ img ] = CTbackprojectionv2( proj, param )
%CTBACKPROJECTION Summary of this function goes here
%   Detailed explanation goes here

img = zeros(param.nPx, param.nPx, param.nPx, 'single');
[x,y] = meshgrid(param.x,param.y);
if param.gpu == 1
        x = gpuArray(x);
    y = gpuArray(y);
end

for i = 1:param.nProj
    disp(i);
    tic
    img = img + backprojection(proj(:,:,i),param,i,x,y);
    t(i) = toc;
end
toc

% figure;
% plot(t);
% title(strcat(num2str(mean(t)),'\pm',num2str(std(t))));

end

function vol = backprojection(proj,param,iview,x,y)

beta = param.theta(iview)/360*2*pi;
vol = gpuArray(zeros(param.nPx,param.nPx,param.nPx,'single'));

if param.off_u<=0
u = single(param.D.*(x*cos(beta)+y*sin(beta)+abs(param.off_u))...
    ./(x*sin(beta)-y*cos(beta)+param.D)-abs(param.off_u));
else
u = single(param.D.*(x*cos(beta)+y*sin(beta)-abs(param.off_u))...
    ./(x*sin(beta)-y*cos(beta)+param.D)+abs(param.off_u));
end

frac = single(param.D./(x*sin(beta)-y*cos(beta)+param.D));

if param.gpu == 1
    proj = gpuArray(single(proj));
    ys = gpuArray(param.y(1)-1);
    u2 = gpuArray(u-param.x(1)+1);
end    

for iz = 1:param.nPx
    if param.gpu == 1
        v = (frac.*(-param.nPx/2-param.opticCentreP(2)+iz-1))+param.off_v;
        v2 = v-ys;
        vol(:,:,iz) = flipud(gather(param.D^2.*frac.*interp2(proj',u2,v2,param.interptype)));
    else
        v = single(frac.*param.z(iz)+param.off_v);
        vol(:,:,iz) = param.D^2.*frac./(param.D^2+v.^2).*interp2(proj',u-param.x(1)+1,v-param.y(1)+1,param.interptype);
    end 
    
    
end
vol(isnan(vol))=0;

return

end

