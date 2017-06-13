function [ optimalShift ] = findCOR( sinogram )
% Written by T.Watson 16/10/2015. Adapted from algorithm presented in
% Reliable Method for calculating centre of rotation in parallel beam
% tomography 2014
k=1;
shiftRange = -50:1:50;
for s = shiftRange

%copy sinogram. sinogram is from 0-180 degrees. Flip and shift and add to
%other half to form full sinogram shifted about horizontal centre.

shiftedSino = circshift(sinogram,[s,0]);

copy = flipud(shiftedSino);

full = [shiftedSino,copy];

% imagesc(full);
% drawnow;

ft = fftshift(fft2(full));

dx = 1;

%max radial extent is half the sinogram x dimension
r = size(full,1)/2;

dtheta = 2*pi/size(full,2);

u = linspace(-1/(2*dx),1/(2*dx),size(full,1));

n = linspace(-1/(2*dtheta),1/(2*dtheta),size(full,2));


[n2,u2] = meshgrid(n,u);

mask = zeros(size(full,1),size(full,2));

mask(or(n2>abs(r.*u2),n2<-abs(r.*u2)))=1;
ft = abs(ft).*mask;

Q(k) = sum(ft(:))/sum(mask(:));

k=k+1;

%% plot options
% figure;
% subplot(1,2,1)
% imagesc(full);
% xlabel('Theta');
% ylabel('x');


% y = r.*u;
% figure;
% plot(y,u);
% hold on
% plot(-y,u);
% hold off
% xlim([min(n),max(n)]);
% ylim([min(u),max(u)]);
% xlabel('n');
% ylabel('u');

% figure;
% subplot(1,2,2)
% imagesc(n,u,abs(ft));
% xlabel('n');
% ylabel('u');
% caxis([0,10e4]);
% drawnow;
% pause(0.1);

end
%figure;
plot(linspace(min(shiftRange),max(shiftRange),length(Q)),Q);
drawnow;
[~,lc] = min(Q);

optimalShift = mean(shiftRange(lc));

end

