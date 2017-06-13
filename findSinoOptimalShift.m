function [ optimalShift ] = findSinoOptimalShift( sinogram )
%Find optimal shift across 360 degrees of projections for sinogram
k=1;
for i=1:3:size(sinogram,2)/2
tic
    sinoPortion = sinogram(:,i:i+size(sinogram,2)/2);
%     imagesc(sinoPortion);
%     drawnow;
    
    optimalShift(k) = findCOR( sinoPortion );
    k=k+1;
    disp(i);
toc
end

figure;
plot(optimalShift);
end

