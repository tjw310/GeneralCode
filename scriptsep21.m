% z1 = zeros(images.nPx); z2 = zeros(images.nPx);
% for i=1:1:400
% n=i;
% s1 = real(testFan(images,n));
% s2 = real(testFan(images2,n));
% z1 = s1+z1; z2 = s2+z2;
% if rem(i,100)==0
%     subplot(2,2,1); imagesc(s1-s2); axis square;
%     subplot(2,2,2); imagesc(z1); axis square;
%     subplot(2,2,3); imagesc(z2); axis square;
%     subplot(2,2,4); imagesc(real(z1-z2)); axis square;
%     drawnow; end
% end
% z1(z1<0) = 0; z2(z2<0)=0; figure; imagesc(z1);
%%
z1 = zeros(images2.nPx); z2 = zeros(images2.nPx);
figure;
for i=1:images2.nProj
s1 = testFan(images2,i); 
s2 = testFanv3(images2,i); 
z1 = s1+z1; z2 = s2+z2;
if rem(i,10)==0
subplot(1,2,1); imagesc(real(z1)); axis square; drawnow; 
subplot(1,2,2); imagesc(real(z2)); axis square; drawnow; end
end