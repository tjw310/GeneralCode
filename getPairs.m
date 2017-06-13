function [x1,x2,z1,z2,y1] = getPairs(images,fixed)

[x,z] = FastPeakFind(images.projections(:,:,1).');

figure; imagesc(images.projections(:,:,1).'); hold on; scatter(x,z,'rx'); hold off;



[x180,z180] = FastPeakFind(images.projections(:,:,201).');

figure; imagesc(images.projections(:,:,201).'); hold on; scatter(x180,z180,'rx'); hold off;

x180 = 2600-x180;

figure;
scatter(x,z,'b'); hold on; scatter(x180,z180,'rx'); hold off; drawnow;

k=1;
loopno = 1;
while ~isempty(x)
    r = sqrt((x(1)-(2600-x180)).^2+(z(1)-z180).^2);
    [mn,lc] = min(r);
    if mn<50
        %z(k,2) = ob1(i).zDepth;
        %z(k,1) = ob2(lc).zDepth;
        %z(k,5) = (ob1(i).zDepth+ob2(lc).zDepth)/2;
        %z(k,3) = ob1(i).centreX;
        %z(k,4) = ob1(i).centreY;
        %z(k,6) = i;
        %z(k,7) = lc;
        
        z1(k) = z(1)-2160/2;
        z2(k) = z180(lc)-2160/2;
        x1(k) = x(1)-2560/2;
        y1(k) = fixed(loopno).zDepth;
        x2(k) = (x180(lc))-2560/2;
        %y2(k) = ob2(lc).zDepth;
        k=k+1;
        
        x180 = x180((1:length(x180))~=lc);
        z180 = z180((1:length(z180))~=lc);
    end
    x = x(2:end);
    z = z(2:end);
    loopno = loopno+1;
end




figure; scatter(x1,z1,'b'); hold on; scatter(-x2,z2,'rx'); hold off;
end