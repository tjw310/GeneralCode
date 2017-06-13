obj = images2;

x = obj.x;
y = obj.y;

[xx,yy] = meshgrid(x,y);


%k = zeros(obj.nPx,obj.nPx);
%m = zeros(obj.nPx,obj.nPx);

ob = ones(obj.nPx);
obj.nProj = 400;
t = obj.theta;
x = rOrbit.*sin((t)/180*pi);
for i=1:obj.nProj
n = round(x(i)/stepsize);
discretemotion(i) = n*stepsize;
end
obj.AoRmotion = discretemotion;
close all; clear k k2 m m2;
for i=1:obj.nProj
    beta = t(i)/360*2*pi;

    pos = obj.AoRcentreX;
    sqrd = obj.opticCentre(1)/obj.SubSampFc;
    OsOa = (pos-obj.opticCentre(1))/obj.SubSampFc;
    
    l = obj.D.*(xx.*cos(beta)+yy.*sin(beta)+OsOa)./(xx.*sin(beta)-yy.*cos(beta)+obj.D)+sqrd;
    l2 = l;
    l(l<min(obj.x))=NaN; l(l>max(obj.x))=NaN;
    l2(and(l2>=min(obj.x),l2<=max(obj.x))) = NaN;
    gamma = atan2(l,obj.D)/pi*180;
    theta = gamma+beta/pi*180;
    taxis = obj.D.*sin(gamma*180/pi);
    p = gamma./180;
   % gamma(p<0) = gamma(p<0)+180;
    %gamma(p>=0)= gamma(p>=0)-180*floor(p(p>=0));
    gamma2 = (atan2(l2,obj.D)+beta)/pi*180;
    p2 = gamma2./180;
   % gamma2(p2<0) = gamma2(p2<0)+180;
    %gamma2(p2>=0)= gamma2(p2>=0)-180*floor(p2(p2>=0));
    k(:,:,i) = gamma;
    k2(:,:,i) = gamma2;
    subplot(1,2,1); imagesc(gamma); axis square; caxis([0,360]);
    
    pos = obj.AoRcentreX+obj.AoRmotion(i);
    OsOa = (pos-obj.opticCentre(1))/obj.SubSampFc; 
    
    l = obj.D.*(xx.*cos(beta)+yy.*sin(beta)+OsOa)./(xx.*sin(beta)-yy.*cos(beta)+obj.D)+sqrd;
    l2 = l;
    l(l<min(obj.x))=NaN; l(l>max(obj.x))=NaN;
    l2(and(l2>=min(obj.x),l2<=max(obj.x))) = NaN;
    gamma = (atan2(l,obj.D)+beta)/pi*180;
    p = gamma./180;
    %gamma(p<1/180) = gamma(p<1/180)+180;
    %gamma(p>=1/180)= gamma(p>=1/180)-180*floor(p(p>=1/180));
    gamma2 = (atan2(l2,obj.D)+beta)/pi*180;
    p2 = gamma2./180;
    %gamma2(p2<0) = gamma2(p2<0)+180;
    %gamma2(p2>=0)= gamma2(p2>=0)-180*floor(p2(p2>=0));
    m(:,:,i) = gamma;
    m2(:,:,i) = gamma2;
    subplot(1,2,2); imagesc(gamma); axis square; 
    caxis([0,360]); 
    drawnow; 
end

%%
figure;
for i=1:512
    r = squeeze(k(i,:,:)).'; subplot(2,2,1); %imagesc(r);
    r2 = squeeze(k2(i,:,:)).'; subplot(2,2,2); r(isnan(r))=0;r2(isnan(r2))=0;
    rSum = (r2+r); %imagesc(rSum);
    [r3,I]=sort(rSum,1,'ascend');%subplot(2,2,3); imagesc(r3);
    [rows,cols] = size(rSum);
    I = I +repmat((0:cols-1),rows,1)*rows;
    r = r(I); subplot(2,2,4); %imagesc(r); drawnow;
    kSort(i,1:size(r,2),1:size(r,1)) = r.';
    r = squeeze(m(i,:,:)).'; subplot(2,2,1); %imagesc(r);
    r2 = squeeze(m2(i,:,:)).'; subplot(2,2,2); r(isnan(r))=0;r2(isnan(r2))=0;
    rSum = (r2+r); %imagesc(rSum);
    [r3,I]=sort(rSum,1,'ascend'); %subplot(2,2,3); imagesc(r3);
    [rows,cols] = size(rSum);
    I = I +repmat((0:cols-1),rows,1)*rows;
    r = r(I); subplot(2,2,4); %imagesc(r); drawnow;
    mSort(i,1:size(r,2),1:size(r,1)) = r.';
end
%%
v = squeeze(k(256,73,:));
v2 = squeeze(k2(256,73,:));
v(isnan(v))=0;
v2(isnan(v2)) = 0;
[~,I]= sort(v+v2,'ascend');
v2 = v2(I);
v = v(I);
w = squeeze(m(256,73,:));
w2 = squeeze(m2(256,73,:));
w(isnan(w))=0;
w2(isnan(w2)) = 0;
[~,I]= sort(w+w2,'ascend');
w2 = w2(I);
w = w(I);
figure; scatter(1:360,v); hold on; scatter(1:360,v2); hold off;
figure; scatter(1:360,w); hold on; scatter(1:360,w2); hold off;
%%
close all
for i=256
    subplot(1,2,1); imagesc(squeeze(k(:,i,:)).'); axis equal tight;
    %subplot(1,2,2); imagesc(squeeze(k(i,:,1:180))+fliplr(flipud(squeeze(k(i,:,181:end))))); colorbar;
    subplot(1,2,2); imagesc(squeeze(m(:,i,:)).'); axis equal tight;
    drawnow;
end

%%
line = floor(squeeze(k(20,20,:)));
[a,b,c ] =unique(line(~isnan(line)));
figure; plot(line); hold on; plot(c); hold off;