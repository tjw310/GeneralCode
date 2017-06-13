function [ img ] = CTbackprojection( proj, param )
%CTBACKPROJECTION Summary of this function goes here
%   Detailed explanation goes here

img = zeros(param.nx, param.ny, param.nz, 'single');

for i = 1:param.nProj
    disp(i);
    img = img + backprojection(proj(:,:,i),param,i);
end


end

function vol = backprojection(proj,param,iview)

angle_rad = param.theta(iview)/360*2*pi;
vol = zeros(param.nx,param.ny,param.nz,'single');

[xx,yy] = meshgrid(param.xs,param.ys);

rx = xx.*cos(angle_rad-pi/2) + yy.*sin(angle_rad-pi/2);
ry = -xx.*sin(angle_rad-pi/2) + yy.*cos(angle_rad-pi/2);

rx = xx.*cos(angle_rad)-yy.*sin(angle_rad);
ry = xx.*cos(angle_rad)+yy.*sin(angle_rad);

pu = single(((rx.*(param.DSD)./(ry + param.DSO))+param.us(1))/(-param.du) + 1);
Ratio = (single(param.DSO.^2./(param.DSO+ry).^2));    
if param.gpu == 1
    pu = gpuArray(single(pu));
    proj = gpuArray(single(proj));
    Ratio = gpuArray(Ratio);
end    

for iz = 1:param.nz   
    if param.gpu == 1
        pv = gpuArray(single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1));
        vol(:,:,iz) = gather(Ratio.*interp2(proj',pu,pv,param.interptype));
    else
        pv = single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1);
        vol(:,:,iz) = (Ratio.*interp2(proj',pu,pv,param.interptype));
    end 
    
    
end

vol(isnan(vol))=0;

return

end

