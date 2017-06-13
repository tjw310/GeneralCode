ob1 = fixed0;
ob2 = fixed180;
%%
clear z1 z2 x1 y1 x2 y2

for i=1:length(ob1)
    x(i) = ob1(i).centreX;
    z(i) = ob1(i).centreY;
end

for i=1:length(ob2)
    x180(i) = ob2(i).centreX;
    z180(i) = ob2(i).centreY;
end

scatter(x,z); hold on; scatter(x180,z180); hold off;

k=1;
loopno=1;

while ~isempty(x)
    r = sqrt((x(1)-(2600-x180)).^2+(z(1)-z180).^2);
    [mn,lc] = min(r);
    if mn<30 && abs((ob1(loopno).zDepth+ob2(lc).zDepth)/2 )<0.02
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
        y1(k) = ob1(loopno).zDepth;
        x2(k) = (x180(lc))-2560/2;
        y2(k) = ob2(lc).zDepth;
        k=k+1;
        
        x180 = x180((1:length(x180))~=lc);
        z180 = z180((1:length(z180))~=lc);
    end
    x = x(2:end);
    z = z(2:end);
    loopno = loopno+1;
end
scatter(z1,(y2+y1)/2); drawnow;

[e,m] = fitLinear(z1,(y2+y1)/2);

[~,fit] = m(e);

hold on; plot(z1,fit); hold off

figure; hold on; scatter(x1,z1,'b'); 
scatter(images.opticCentre(1),images.opticCentre(2),'k');
scatter(x1,z1,[],y1,'filled'); scatter(-x2+50,z2,'r');  scatter(-x2+50,z2,[],y2,'filled');hold off;
colormap(jet);

% [~,I] = sort(z(:,6),'ascend');
% 
% for i=1:size(z,2)
%    z2(:,i) = z(I,i);
% end
% 
% figure; scatter(z2(:,6),z2(:,4));
%%
% ob1 = mapped0;
% ob2 = mapped180;
% k=1;
% for i=1:length(ob1)
%     r=[];
%     for j=1:length(ob2)
%         r(j) = sqrt((ob1(i).statX-2560+ob2(j).statX).^2+(ob1(i).statY-ob2(j).statY).^2);
%     end
%     [mn,lc] = min(r);
%     if mn<50
%         ob1(1,i)
%         pairs(k,1) = ob1(1,i);
%         pairs(k,2) = ob2(1,lc);
%         %zDepths(k,1) = ob1(1,i).zDepth;
%         %zDepths(k,2) = ob2(1,lc).zDepth;
%         %zDepths(k,3) = (ob1(1,i).zDepth+ob2(1,lc).zDepth)/2;
%         k=k+1;
%     end
% end
%         