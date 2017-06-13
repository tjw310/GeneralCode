function reconConvent(obj,finalsize,imshow,varargin)
%imshow==1 to display image
if nargin>3
    if nargin==5
        stIm = varargin{1};
        fnIm = varargin{2};
    elseif nargin==4
        stIm = varargin{1};
        fnIm = varargin{1};
    else
        error('Enter either slice no, or stIm and fnIm');
    end
else
    stIm = 1;
    fnIm = size(obj.projections,2);
end
        

if obj.AoRangle~=0
    for i=1:obj.nProj
        p = obj.projections(:,:,i);
        p = imrotate(p,obj.AoRangle,'crop');
        obj.projections(:,:,i) = p;
    end
end

ctf = round((size(obj.projections,1)-finalsize)/2);

for i=stIm:fnIm
    sino = gpuArray(squeeze(obj.projections(ctf+1-obj.AoRcentreX:size(obj.projections,1)-ctf-obj.AoRcentreX,i,:)));
    sino(isnan(sino))=0;
    slice = iradon(sino,linspace(0,360-360/400,400),'linear','Ram-lak',1,size(sino,1));
    slice = slice+10000;
    %slice = slice*1000;
    if imshow==1
        imagesc(slice); drawnow;
    end
    disp((i-stIm)/(fnIm-stIm)*100)
    slice = gather(uint16(slice));
    if i==stIm
    name = 'reconSlices';
    mkdir(obj.path,name);
    p = fullfile(obj.path,name);
    end
    imwrite(slice,strcat(fullfile(p,num2str(i,'%05d')),'.tiff'));
end

end