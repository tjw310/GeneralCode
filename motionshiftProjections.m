function out = motionshiftProjections(obj,stageamplitude)

xmotion = getXmotion(obj,stageamplitude)/obj.pxSz*obj.mtz0;

[x,y]=meshgrid(1:size(obj.projections,1));

for i=1:size(obj.projections,3)    
    out(:,:,i) = interp2(obj.projections(:,:,i),y,x+xmotion(i)); 
    disp(i/size(obj.projections,3));
end
end