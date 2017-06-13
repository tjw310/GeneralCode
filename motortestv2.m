function [path,progX,progZ,xOut,zOut,estimates] = motortestv2(varargin)

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
    
    [x,z] = FastPeakFind(im,1);
    
    imagesc(im); hold on; scatter(x,z,'rx'); 
    if i==1 
        xInit=[]; zInit=[];
        c = input('enter coord (press blank if fin): ');
        while ~isempty(c)
            xInit = horzcat(xInit,c(1));
            zInit = horzcat(zInit,c(2));
            c = input('enter coord (press lbank if finished): ');
        end
        for j=1:length(xInit)
            r = sqrt((x-xInit(j)).^2+(z-zInit(j)).^2);
            [~,lc] = min(r);
            xOut(k,j) = x(lc);
            zOut(k,j) = z(lc);
        end
    else
        xPred = repmat(-(progX(k)-progX(k-1))./pxSz*mgGuess,1,length(xInit))+xOut(k-1,:);
        zPred = repmat(-(progZ(k)-progZ(k-1))./pxSz*mgGuess,1,length(xInit))+zOut(k-1,:);
        for j=1:length(xInit)
            r = sqrt((x-xPred(j)).^2+(z-zPred(j)).^2);
            [~,lc] = min(r);
            xOut(k,j) = x(lc);
            zOut(k,j) = z(lc);
        end
    end   
    scatter(xOut(k,:),zOut(k,:),'g'); hold off; drawnow;
    k=k+1;
end

estimates = fitAugust31(xOut,zOut,progX);  
title(sprintf('mag=%.3f',abs(6.5e-3*estimates(1))));
savefig(fullfile(path,sprintf('mag=%.3f.fig',abs(6.5e-3*estimates(1)))));

end