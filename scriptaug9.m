t = 360-images.theta;
amp = 0.1;
motx = amp*sin(t/180*pi);
moty = amp*cos(t/180*pi);

%motx = ones(1,length(t)).*1;
close all
testTraceDif(0,1,motx,moty,images);
testTraceDif(5,0,motx,moty,images);
testTraceDif(5,0.001,motx,moty,images);