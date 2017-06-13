function [xmotion] = getXmotion(obj,stageamplitude)

y = obj.piezomotion;
e = fitSinusoid(y);
xmotion = stageamplitude.*sin(obj.theta/180*pi+e(2)-pi/2);

end