function [path,progX,progZ,xOut,zOut,estimates] = stagetestmotion(varargin)

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
for i=1:nProj
    imInfo = strsplit(im_dir(i).name,{'x=','z=','y=','.t'});
    n(i) = str2num(cell2mat(imInfo(1)));
end

[~,sortI] = sort(n,'ascend');
figure;
mgGuess = 20.27;
pxSz = 6.5e-3;
k=1;
for i=sortI
    disp(k/nProj*100);
    im = single(imread(strcat(path,'\',char(cellstr(im_dir(i).name)))));
    imInfo = strsplit(im_dir(i).name,{'x=','z=','y=','.t'})
    
    progX(k) = str2num(cell2mat(imInfo(3)));
    progZ(k) = str2num(cell2mat(imInfo(2)));
    
    if i==1
        yLoc = str2num(cell2mat(imInfo(4)));
    end
    
    [x,z] = FastPeakFind(im,1);
    
    
    imagesc(im); hold on; scatter(x,z,'rx'); 
    if i==1
        xInit = input('enter x-loc: ');
        zInit = input('enter z-loc: ');
        r = sqrt((x-xInit).^2+(z-zInit).^2);
        [~,lc] = min(r);
        xOut(k) = x(lc);
        zOut(k) = z(lc);
    else
        xPred = -(progX(k)-progX(k-1))./pxSz*mgGuess+xOut(k-1);
        zPred = -(progZ(k)-progZ(k-1))./pxSz*mgGuess+zOut(k-1);
        r = sqrt((x-xPred).^2+(z-zPred).^2);
        [~,lc] = min(r);
        xOut(k) = x(lc);
        zOut(k) = z(lc);
    end   
    scatter(xOut(k),zOut(k),'g'); hold off; drawnow;
    k=k+1;
end

estimates = fitStageMotionToProgValues(xOut,zOut,progX);  
mag = estimates(1)*6.5e-3;
title(sprintf('x=%d,z=%d,y=%d,mag is %d',xInit,zInit,yLoc,mag));
savefig(fullfile(path,sprintf('x=%d,z=%d,y=%.3f,mag is %.2f.fig',xInit,zInit,yLoc,mag)));

end