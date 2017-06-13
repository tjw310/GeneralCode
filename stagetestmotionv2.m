function [path,mipRed,mipGreen,yRed,yGreen,tform,stpR,stpG,magR,magG] = stagetestmotionv2(varargin)

if nargin>0
    path = varargin{1};
else
    path = uigetdir();
end

im_dir = dir(fullfile(path,'*.tif'));
if isempty(im_dir)
    im_dir = dir(fullfile(path,'*.tiff'));
end

nProj = length(im_dir);
figure;
mgGuess = 50;
pxSz = 6.5e-3;
k=1;

for i=1:2:nProj
    disp(k/nProj*200);
    
    info1 = strsplit(im_dir(i).name,{'x=','z=','y=','ETL=','mA','.t'});
    info2 = strsplit(im_dir(i+1).name,{'x=','z=','y=','ETL=','mA','.t'});
    
    % The 'moving' picture is chosen to be the one with the lower y= /
    % ETL= number. NOTE this will give opposite magnifications for ETL and
    % piezo.
    if str2num(cell2mat(info2(4)))>str2num(cell2mat(info1(4)))
        imRed = single(imread(strcat(path,'\',char(cellstr(im_dir(i).name)))));
        imInfoRed = strsplit(im_dir(i).name,{'x=','z=','y=','ETL=','mA','.t'});
        progX(k) = str2num(cell2mat(imInfoRed(2)));
        progZ(k) = str2num(cell2mat(imInfoRed(3)));
        yRed(k) = str2num(cell2mat(imInfoRed(4)));

        imGreen = single(imread(strcat(path,'\',char(cellstr(im_dir(i+1).name)))));
        imInfoGreen = strsplit(im_dir(i+1).name,{'x=','z=','y=','ETL=','mA','.t'});
        progXgreen(k) = str2num(cell2mat(imInfoGreen(2)));
        progZgreen(k) = str2num(cell2mat(imInfoGreen(3)));
        yGreen(k) = str2num(cell2mat(imInfoGreen(4)));   
    else
        imRed = single(imread(strcat(path,'\',char(cellstr(im_dir(i+1).name)))));
        imInfoRed = strsplit(im_dir(i+1).name,{'x=','z=','y=','ETL=','mA','.t'});
        progX(k) = str2num(cell2mat(imInfoRed(2)));
        progZ(k) = str2num(cell2mat(imInfoRed(3)));
        yRed(k) = str2num(cell2mat(imInfoRed(4)));

        imGreen = single(imread(strcat(path,'\',char(cellstr(im_dir(i).name)))));
        imInfoGreen = strsplit(im_dir(i).name,{'x=','z=','y=','ETL=','mA','.t'});
        progXgreen(k) = str2num(cell2mat(imInfoGreen(2)));
        progZgreen(k) = str2num(cell2mat(imInfoGreen(3)));
        yGreen(k) = str2num(cell2mat(imInfoGreen(4)));
    end
    
    if i==1; thres = mean(imRed(:)); end
    
    [xRed,zRed] = FastPeakFindv2(imRed,1,thres);
    [xGreen,zGreen] = FastPeakFindv2(imGreen,1,thres); 
    if i==1
        imagesc(imRed); hold on; scatter(xRed,zRed,'rx'); hold off;
        xInit = input('enter Red x-loc: ');
        zInit = input('enter Red z-loc: ');
        r = sqrt((xRed-xInit).^2+(zRed-zInit).^2);
        [~,lc] = min(r);
        xOut(k) = xRed(lc);
        zOut(k) = zRed(lc);
        
        imagesc(imGreen); hold on; scatter(xGreen,zGreen,'rx'); hold off;
        xInit = input('enter Green x-loc: ');
        zInit = input('enter Green z-loc: ');
        r = sqrt((xGreen-xInit).^2+(zGreen-zInit).^2);
        [~,lc] = min(r);
            if ~isempty(lc)
                xOutG(k) = xGreen(lc);
                zOutG(k) = zGreen(lc);

                av = mean(imRed(:));
                mipRed = ones(size(imRed)).*av;
                mipGreen = ones(size(imRed)).*mean(imGreen(:));
                idxRedz = zOut(k)-60:zOut(k)+60; idxRedz(or(idxRedz<1,idxRedz>size(imRed,1))) = NaN;
                idxRedx = xOut(k)-60:xOut(k)+60; idxRedx(or(idxRedx<1,idxRedx>size(imRed,2))) = NaN;
                idxRedz = idxRedz(~isnan(idxRedz)); idxRedx = idxRedx(~isnan(idxRedx));
                mipRed(idxRedz,idxRedx) = imRed(idxRedz,idxRedx)./max(max(imRed(idxRedz,idxRedx)));
                idxGreenz = zOutG(k)-60:zOutG(k)+60; idxGreenz(or(idxGreenz<1,idxGreenz>size(imGreen,1))) = NaN;
                idxGreenx = xOutG(k)-60:xOutG(k)+60; idxGreenx(or(idxGreenx<1,idxGreenx>size(imGreen,2))) = NaN;
                idxGreenz = idxGreenz(~isnan(idxGreenz)); idxGreenx = idxGreenx(~isnan(idxGreenx));
                mipGreen(idxGreenz,idxGreenx) = imGreen(idxGreenz,idxGreenx)./max(max(imGreen(idxGreenz,idxGreenx)));
                k=k+1;
            end
    else
        %if k==2; dpx = abs(-(progX(k)-progX(k-1))./pxSz*mgGuess); end
        xPred = -(progX(k)-progX(k-1))./pxSz*mgGuess+xOut(k-1);
        zPred = -(progZ(k)-progZ(k-1))./pxSz*mgGuess+zOut(k-1);
        r = sqrt((xRed-xPred).^2+(zRed-zPred).^2);
        [~,lc] = min(r);
        xOut(k) = xRed(lc);
        zOut(k) = zRed(lc);
        
        xPred = -(progXgreen(k)-progXgreen(k-1))./pxSz*mgGuess+xOutG(k-1);
        zPred = -(progZgreen(k)-progZgreen(k-1))./pxSz*mgGuess+zOutG(k-1);
        r = sqrt((xGreen-xPred).^2+(zGreen-zPred).^2);
        [~,lc] = min(r);
            if ~isempty(lc)
            xOutG(k) = xGreen(lc);
            zOutG(k) = zGreen(lc);

            %subplot(1,2,2); imagesc(imGreen); hold on; scatter(xPred,zPred,'w'); hold off;

            idxRedz = zOut(k)-60:zOut(k)+60; idxRedz(or(idxRedz<1,idxRedz>size(imRed,1))) = NaN;
            idxRedx = xOut(k)-60:xOut(k)+60; idxRedx(or(idxRedx<1,idxRedx>size(imRed,2))) = NaN;
            idxRedz = idxRedz(~isnan(idxRedz)); idxRedx = idxRedx(~isnan(idxRedx));
            mipRed(idxRedz,idxRedx) = imRed(idxRedz,idxRedx);
            idxGreenz = zOutG(k)-60:zOutG(k)+60; idxGreenz(or(idxGreenz<1,idxGreenz>size(imGreen,1))) = NaN;
            idxGreenx = xOutG(k)-60:xOutG(k)+60; idxGreenx(or(idxGreenx<1,idxGreenx>size(imGreen,2))) = NaN;
            idxGreenz = idxGreenz(~isnan(idxGreenz)); idxGreenx = idxGreenx(~isnan(idxGreenx));
            mipGreen(idxGreenz,idxGreenx) = imGreen(idxGreenz,idxGreenx);

            k=k+1;
            end
    end
%     subplot(1,2,1); imagesc(mipRed); axis square; hold on;
%     scatter(xOut(k),zOut(k),'g'); hold off;
%     subplot(1,2,2); imagesc(mipGreen); axis square;
%     drawnow;
  %  k=k+1;
end

stepSize = progZ(1:end-1)-progZ(2:end); stepSize(stepSize==0)=NaN;
stepSize = nanmean(stepSize);

[tform,gradError,~,~,~,~,~,~,stpR,stpG] = guessInitialTransform(mipRed,mipGreen);
movReg = imwarp(mipRed,tform,'outputview',imref2d(size(mipGreen)));
figure; subplot(1,2,1); imshowpair(mipGreen,mipRed,'scaling','joint');
title(sprintf('deltaM=%.4f',tform.T(1,1)));
subplot(1,2,2); imshowpair(mipGreen,movReg,'scaling','joint'); drawnow;
title(sprintf('errdeltaM=%.6f',gradError));
savefig(fullfile(path,sprintf('deltaM=%.4f,err=%.6f.fig',tform.T(1,1),gradError)));

magR = stpR./(stepSize/pxSz);
magG = stpG./(stepSize/pxSz);

figure; plot(magR); hold on; plot(magG); hold off;
legend('Red Mag','Green Mag');
xlabel('Step Number'); ylabel('Magnification');
title(sprintf('Step Magnification'));
savefig(fullfile(path,sprintf('Step Magnification.fig')));

end