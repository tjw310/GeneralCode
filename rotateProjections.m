function projOut = rotateProjections(obj)

[X,Z]= meshgrid(obj.x,obj.z);
if obj.AoRangle~=0
    X = X*cos(obj.AoRangle/180*pi)-Z*sin(obj.AoRangle/180*pi)-X(1,1);
    Z = Z*cos(obj.AoRangle/180*pi)+X*sin(obj.AoRangle/180*pi)-Z(1,1);
end
projOut = zeros(obj.nPx,obj.nPx,obj.nProj,'single');
    for i=1:obj.nProj
        disp(i/obj.nProj*100)

        single_scaled_proj = single(interp2(obj.projections(:,:,i).',X,Z,obj.interptype));
        %subplot(1,2,1); imagesc(obj.projections(:,:,i).'); subplot(1,2,2); imagesc(single_scaled_proj); drawnow;
        projOut(:,:,i) = single_scaled_proj.';  
    end        
end