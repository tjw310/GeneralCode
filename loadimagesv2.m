
function [images,path] = loadimagesv2()
% reads in images from folder. Must be .tiff or .tif
path = uigetdir('Z:\','Select Folder of Images');

im_dir = dir(fullfile(path,'*.tiff'));
if size(im_dir,1)==0
    im_dir = dir(fullfile(path,'*.tif'));
end
if size(im_dir)==0
    error('No .tiff images in selected folder');
end

im_size = size(imread(strcat(path,'\',char(cellstr(im_dir(1).name)))));

images = single(zeros(max(im_size),min(im_size),length(im_dir)));

for i=1:length(im_dir)
    
    if im_size(1,1)>im_size(1,2)
        images(1:max(im_size),1:min(im_size),i) = (single(imread(strcat(path,'\',char(cellstr(im_dir(i).name))))));
    else
        images(1:max(im_size),1:min(im_size),i) = transpose(single(imread(strcat(path,'\',char(cellstr(im_dir(i).name))))));
    end
    disp(i/length(im_dir)*100);
end

end
    




