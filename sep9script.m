gamma = images2.AoRangle/180*pi;
op = images.RopticCentre;
t1 = t;

%opR(1) = op(1).*cos(gamma)-op(2).*sin(gamma);
%opR(2) = op(2).*cos(gamma)+op(1).*sin(gamma);

%t2(1,:,:) = t(1,:,:).*cos(gamma)-t(3,:,:).*sin(gamma);
%t2(2,:,:) = t(2,:,:);
%t2(3,:,:) = t(3,:,:).*cos(gamma)+t(1,:,:).*sin(gamma);

%dmdy = images.dmdy;

%[t1(1,:,:),t1(3,:,:)] = scaleCoords(t1(1,:,:),t1(2,:,:),t1(3,:,:),op,dmdy);
%[t2(1,:,:),t2(3,:,:)] = scaleCoords(t2(1,:,:),t2(2,:,:),t2(3,:,:),opR,dmdy);
        
t3(1,:,:) = t1(1,:,:).*cos(gamma)-t1(3,:,:).*sin(gamma);
t3(2,:,:) = t1(2,:,:);
t3(3,:,:) = t1(3,:,:).*cos(gamma)+t1(1,:,:).*sin(gamma);

figure; subplot(2,2,1); plot(squeeze(t2(1,2,:))); hold on; plot(squeeze(t3(1,2,:))); hold off;
subplot(2,2,2); plot(squeeze(t2(2,2,:))); hold on; plot(squeeze(t3(2,2,:))); hold off;
subplot(2,2,3); plot(squeeze(t2(3,2,:))); hold on; plot(squeeze(t3(3,2,:))); hold off;