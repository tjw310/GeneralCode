function [xP,yP,xE,yE]=findScanRange(piezoscanpath,etlscanpath)

im_dir1 = dir(fullfile(piezoscanpath,'*.tif'));
if isempty(im_dir1)
    im_dir1 = dir(fullfile(piezoscanpath,'*.tiff'));
end

no1 = length(im_dir1);

[~,fixedIstr,~] = fileparts(piezoscanpath);
sepStr = strsplit(fixedIstr,{'mA'});
fixedI = sepStr{1};

im_dir2 = dir(fullfile(etlscanpath,'*.tif'));
if isempty(im_dir2)
    im_dir2 = dir(fullfile(etlscanpath,'*.tiff'));
end

[~,fixedZstr,~] = fileparts(etlscanpath);
sepStr = strsplit(fixedZstr,{'mA'});
fixedZ = sepStr{1};

no2 = length(im_dir2);

for i=1:no1
    im = single(imread(strcat(piezoscanpath,'\',char(cellstr(im_dir1(i).name))))).';
    [x,y] = FastPeakFind(im);
    xP(i,1:length(x)) = x;
    yP(i,1:length(x)) = y;
end

for i=1:no1
    im = single(imread(strcat(etlscanpath,'\',char(cellstr(im_dir2(i).name))))).';
    [x,y] = FastPeakFind(im);
    xE(i,1:length(x)) = x;
    yE(i,1:length(x)) = y;
end

end