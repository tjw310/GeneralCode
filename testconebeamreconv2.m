%Test cone beam reconstruction
%% Recon case 1 - Analytic reconstruction: filtered backprojection
% filter='ram-lak','shepp-logan','cosine', 'hamming', 'hann' : (ramp + additional filter)

param = setParametersv3();

proj_filtered = filteringv3(projections,param);
Reconimg = CTbackprojectionv3(proj_filtered, param);

% for i=1:param.nz
%     figure(2); imagesc(max(Reconimg(:,:,i),0)); axis off; axis equal; colormap gray; colorbar;
%     title(num2str(i));
%     pause(0.01);
% end

%%
figure;
subplot(2,2,1)
imagesc(param.rx,param.rz,squeeze(max(Reconimg,[],1)).'); axis square; xlabel('x'); ylabel('z'); hold on;
set(gca,'ydir','normal');
scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); hold off;
subplot(2,2,2)
imagesc(param.ry,param.rz,squeeze(max(Reconimg,[],2)).'); axis square; xlabel('y'); ylabel('z');hold on;
set(gca,'ydir','normal');
scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); hold off;
subplot(2,2,3)
imagesc(param.rx,param.ry,squeeze(max(Reconimg,[],3))); axis square; xlabel('x'); ylabel('y');hold on;
set(gca,'ydir','normal');
scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); hold off;


% %%
% clear slice
% for i=1:size(projections,2);
%     sinogram = squeeze(projections(:,i,:));
%     slice(:,:,i) = iradon(sinogram,param.theta,'ram-lak',1,size(sinogram,1));
% end
% 
% figure;
% subplot(2,2,1)
% imagesc(param.rx,param.rz,(squeeze(max(slice,[],1))).'); axis square; set(gca,'ydir','normal'); hold on;
% scatter(coordsInit(1,:),coordsInit(3,:),'rx'); xlabel('x');ylabel('z'); hold off;
% subplot(2,2,2)
% imagesc(param.ry,param.rz,(squeeze(max(slice,[],2))).');axis square; set(gca,'ydir','normal'); hold on;
% scatter(coordsInit(2,:),coordsInit(3,:),'rx'); xlabel('y');ylabel('z'); hold off;
% subplot(2,2,3)
% imagesc(param.rx,param.ry,squeeze(max(slice,[],3)));axis square; set(gca,'ydir','normal'); hold on;
% scatter(coordsInit(1,:),coordsInit(2,:),'rx'); xlabel('x');ylabel('y'); hold off;
