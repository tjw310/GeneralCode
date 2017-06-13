function projOut = scaleProjections(obj)
[X,Z]= meshgrid(obj.x,obj.z);
if obj.AoRangle~=0
    X = X*cos(obj.AoRangle/180*pi)-Z*sin(obj.AoRangle/180*pi);
    Z = Z*cos(obj.AoRangle/180*pi)+X*sin(obj.AoRangle/180*pi);
    op(1) = obj.opticCentre(1)*cos(obj.AoRangle/180*pi)-obj.opticCentre(2)*sin(obj.AoRangle/180*pi);
    op(2) = obj.opticCentre(2)*cos(obj.AoRangle/180*pi)+obj.opticCentre(1)*sin(obj.AoRangle/180*pi);
else
    op(1) = obj.opticCentre(1);
    op(2) = obj.opticCentre(2);
end
scFc = 1./(1-obj.dmtdy*(1.3330.*obj.piezomotion));
projOut = zeros(obj.nPx,obj.nPx,obj.nProj,'single');
    for i=1:obj.nProj
        disp(i/obj.nProj*100)
        sc = scFc(i);
        X2 = (X-op(1)/obj.SubSampFc)*sc+op(1)/obj.SubSampFc-X(1,1);
        Z2 = (Z-op(2)/obj.SubSampFc)*sc+op(2)/obj.SubSampFc-Z(1,1);

        single_scaled_proj = single(interp2(obj.projections(:,:,i).',X2,Z2,obj.interptype));
        %subplot(1,2,1); imagesc(obj.projections(:,:,i).'); subplot(1,2,2); imagesc(single_scaled_proj); drawnow;
        projOut(:,:,i) = single_scaled_proj.';  
    end        
end