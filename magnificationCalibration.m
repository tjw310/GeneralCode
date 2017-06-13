% script to read in images from labview acquision
% RunTelecentricityCalibration. Access saved folder. Then run algorithm to
% map magnificaiton as a function of piezo offset and ETL scan range
% position.

function sz = magnificationCalibration()

path = uigetdir('Z:\','Select Folder of Images');
fixedpath = fullfile(path,'145mA');
dynamicpath = fullfile(path,'0-290mA');

sz = getsize(fixedpath);



%% get size images folder function
function sz = getsize(path)

im_dir = dir(fullfile(path,'*.tiff'));

if size(im_dir,1)==0
    im_dir = dir(fullfile(path,'*.tif'));
end
if size(im_dir)==0
    error('No .tiff images in selected folder');
end

sz = size(imread(strcat(path,'\',char(cellstr(im_dir(1).name)))));

sz(1,3) = size(im_dir);


end

end
    
