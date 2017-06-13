classdef magnificationCalibrationClassArchive < handle
    %Class and methods for the magnification calibration
    properties
        fixedpath
        dynamicpath
        imagenostart = []
        imagenofinish = []
        mxIm = []
    end
    
    methods
        %% constructors and get/set
        %contructor
        function obj = magnificationCalibrationClass()
            path = uigetdir('Z:\','Select Folder of Images');
            obj.fixedpath = fullfile(path,'145mA');
            obj.dynamicpath = fullfile(path,'0-290mA');
        end
        
        %change path constructor
        function changePath(obj)
            path = uigetdir('Z:\','Select Folder of Images');
            obj.fixedpath = fullfile(path,'145mA');
            obj.dynamicpath = fullfile(path,'0-290mA');
        end
        
        %% public functions
        %get size of images in path, variable is either fixed path or dynamic path
        function sz = sz(obj,type)
            path = obj.testpath(type);
            im_dir = dir(fullfile(path,'*.tiff'));
            if size(im_dir,1)==0
                im_dir = dir(fullfile(path,'*.tif'));
            end
            if size(im_dir)==0
                error('No .tiff images in selected folder');
            end
            sz = size(imread(strcat(path,'\',char(cellstr(im_dir(1).name)))));
            sz(1,3) = length(im_dir);
        end
    
        %load image
        function [image,zDepth] = loadimage(obj,type,imageno)           
            path = obj.testpath(type);
            im_dir = dir(fullfile(path,'*.tiff'));
            imInfo = strsplit(im_dir(imageno).name,{'=','.tiff'});
            zDepth =str2num(imInfo{2});
            image = double(imread(strcat(path,'\',char(cellstr(im_dir(imageno).name)))));
        end
        
        %get beads from image 
        function beads = getBeads(obj,type,imageno)
            [image,zDepth] = loadimage(obj,type,imageno);
            [x,y] = magnificationCalibrationClass.findPeaks(image,[]);
            uPoints(:,1) = y;
            uPoints(:,2) = x;
            beads =obj.testPoints(uPoints,image,zDepth);
            figure;
            magnificationCalibrationClass.scatterPoints(beads,image);
        end
        
        %get beads from multiple images
        function beadsOut = runOverMultipleImages(obj,type,imagenostart,imagenofinish)
            obj.imagenostart = imagenostart;
            obj.imagenofinish = imagenofinish;
            beads = [];
            obj.mxIm=[];
            for i=1:(imagenofinish-imagenostart)
                nBeads = getBeads(obj,type,imagenostart+i-1);
                for j=1:size(nBeads,2)
                    if max(nBeads(j).centIntensity)>0.1*max(obj.mxIm(:))
                    if ~isempty(beads)
                    for k=1:size(beads,2)
                    %find nearest old bead to new bead, and if within 10
                        %pixels add information         
                        dif(k) = sqrt((beads(k).centreX(end)-nBeads(j).centreX).^2+(beads(k).centreY(end)-nBeads(j).centreY).^2);
                    end
                    [vl,lc]=min(dif);
                    if vl<50
                        if and(beads(lc).zDepth~=nBeads(j).zDepth,abs(beads(lc).zDepth(end)-nBeads(j).zDepth)<0.075)
                            addPoint(beads(lc),nBeads(j));
                        end
                    else
                        beads = horzcat(beads,nBeads(j));
                    end
                    
                    else
                        beads = horzcat(beads,nBeads(j));
                    end
                    end
                end
                      
                disp(i/(imagenofinish-imagenostart)*100);
                
                ims(:,:,1) = obj.mxIm;
                newIm = loadimage(obj,type,imagenostart+i-1);
                ims(:,:,2) = newIm;
                obj.mxIm = max(ims,[],3);
            end
            
            k=1;
            for j=1:size(beads,2)
                if length(beads(j).centreX)>5
                    beadsOut(k) = beads(j);
                    k=k+1;
                end
            end
            
            obj.scatterPoints(beadsOut,0);
        end
        
        %function to plot bead centre intensities vs zDepth
        function plotDepthIntensity(obj,beads)
            for i=1:size(beads,2)
                figure;
                scatter(beads(i).zDepth,beads(i).centIntensity);
                %scatter3(beads(i).centreX,beads(i).centreY,beads(i).zDepth,[],beads(i).centIntensity,'filled');
            end
        end
        
        %function to scatter bead centres against maximum of image
        function scatterAllPoints(obj,allbeads,thresholdfraction)
            figure;
            imagesc(obj.mxIm);
            colormap(parula);
            
            hold on
            x=[];
            y=[];
            I=[];
            for i=1:size(allbeads,2)
                I = horzcat(I,allbeads(i).centIntensity);
                x = horzcat(x,allbeads(i).centreX);
                y = horzcat(y,allbeads(i).centreY);
            end
            threshold = thresholdfraction.*max(I);

            I = I(I>threshold);

            x = x(I>threshold);
            y = y(I>threshold);
            scatter(x,y,'r');
            hold off
        end          
        
        %function to find maximum in bead traces and output its location,
        %maximum peak intensity and z depth
        function beadsOut = findBeadZlocations(obj,beads)
            for i=1:size(beads,2)
                [I,lc]=max(beads(i).centIntensity);
                x = beads(i).centreX(lc);
                y = beads(i).centreY(lc);
                z = beads(i).zDepth(lc);
                beadsOut(i) = beadClass(I,x,y,z);
            end
        end
        
        function beadsOut = getFixedBeads(obj)
            path = obj.testpath('fixed');
            im_dir = dir(fullfile(path,'*.tiff'));
            allbeads = obj.runOverMultipleImages('fixed',1,size(im_dir,1));
            beadsOut = obj.findBeadZlocations(allbeads);
        end
        
    end
    %% static methods
    methods (Static)
        %find points via FastPeakFind.m. This uses a local thresholding
        %method. Adapted from mathworks.
        function [x,y] = findPeaks(image,varargin)
            [x,y] = FastPeakFind(image);
            if nargin >1
                figure;
                imagesc(image);
                hold on
                scatter(x,y,'r');
                hold off
            end
        end
        
        %find points above threshold
        function [r,c] = findpoints(image,threshold)
            [r,c] = find(image>threshold);
        end
        
        
        %This gets unique points, defined as points further than 25 pixels from any other point.
        function uPoints = uniquePoints(r,c)
            k=1;
            uPoints=[];
            while and(size(r,1)>0,size(c,1)>0)
                difR = abs(repmat(r(1),size(r,1),1) - r);
                difC = abs(repmat(c(1),size(c,1),1) - c);
                rOut = r(and(difR<25,difC<25));
                cOut = c(and(difR<25,difC<25));
                uPoints(k,1:2) = round(mean([rOut,cOut],1));
                k=k+1;
                r = r(or(difR>=25,difC>=25));
                c = c(or(difR>=25,difC>=25));
            end
        end
        
        %Get unique points via using muliple thresholds and a cluster
        %algorithm
        function uPoints = uniquePointsv2(image,threshold,noThresLevels,areaWidth)
            disp('function.uniquePointsv2');
            thresLevels = linspace(max(image(:)),threshold,noThresLevels);
            R = [];
            C = [];
            figure;
            imagesc(image);
            for i=1:noThresLevels
                threshold = thresLevels(i);
                [r,c]=find(image>threshold);
                if ~isempty(r)
                    %This is the cluster functions, that will try and
                    %minimise the number of clusters as well as there
                    %overlap
                    [nP,rowCentres,colCentres]=clusterPoints(r',c',areaWidth);                    
                    %Ready for mask to hide already found clusters.
                    [x2,y2] = meshgrid(linspace(-1,1,2*areaWidth));
                    r2 = sqrt(x2.^2+y2.^2);
    
                    for k=1:nP
                        mask = ones(2*areaWidth,2*areaWidth);
                        mask(r2<=1) = 0;
                        a = round(rowCentres(k)-areaWidth:rowCentres(k)+areaWidth-1);
                        b = round(colCentres(k)-areaWidth:colCentres(k)+areaWidth-1);
                        mask = mask(~or(a<1,a>2160),~or(b<1,b>2560));
                        a = a(~or(a<1,a>2160));
                        b = b(~or(b<1,b>2560));
                        image(a,b)=image(a,b).*mask;  
                        %drawing
                        %imagesc(image);
                        %hold on
                        %scatter(colCentres,rowCentres,'r');
                        %rectangle('Position',[colCentres(k)-areaWidth/2,rowCentres(k)-areaWidth/2,areaWidth,areaWidth]);
                        %hold off
                    end
                R = horzcat(R,rowCentres);
                C = horzcat(C,colCentres);
                end
            end
            uPoints = [R.',C.'];
            hold on
            scatter(C,R,'r');
            hold off
            
        end

        
        %contour unique points in 50^2 pixel area, to test whether they are
        %beads or not and returns the beads as bead objects, with centres
        %and max intensities
        function beads = testPoints(uPoints,image,zDepth)
            disp('function.testPoints');
            areaWidth = 50;
            beads = [];
            for i=1:size(uPoints,1)
                rows = round(uPoints(i,1)-areaWidth/2:uPoints(i,1)+(areaWidth-1)/2);
                cols = round(uPoints(i,2)-areaWidth/2:uPoints(i,2)+(areaWidth-1)/2);
                rows = rows(and(rows>=1,rows<=2160));
                cols = cols(and(cols>=1,cols<=2560));
                imArea = medfilt2(image(rows,cols),[7 7]);
                %imagesc(imArea);
                %pause(0.3);
                %hold on
                %figure;
                %contour(imArea,10);
                %hold off
                %drawnow;
                c = contourc(imArea,10);
                disp(i);
                if ~isempty(c)
                contours = magnificationCalibrationClass.getcontourlines(c);
                nBeads = magnificationCalibrationClass.sortContour(image,contours,zDepth,min(rows),min(cols));
                if size(nBeads,2)==1
                    if isempty(nBeads.centreX)
                    else
                        beads = horzcat(beads,nBeads);
                    end
                else
                    beads = horzcat(beads,nBeads);
                end
                end
            end 
        end
        
        % function to  get mean and std of contour radii.
        % Then pass through to function testing if there are beads in that
        % portion and how many beads.
        function beads = sortContour(image,contours,zDepth,rowsStart,colsStart)
            k=1;
            params=[];
            beads = [];
            for i=1:size(contours,2)
                   x = contours(1,i).x;
                   y = contours(1,i).y;
                   if length(x)>3
                       r = sqrt((x-mean(x)).^2+(y-mean(y)).^2);
                       params(1,k) = mean(r);
                       params(2,k) = std(r);
                       y2 = round(mean(y))+rowsStart;
                       x2 = round(mean(x))+colsStart;
                       x2(x2<1)=1;
                       y2(y2<1)=1;
                       x2(x2>2560)=2560;
                       y2(y2>2160)=2160;
                       params(3,k) = image(y2,x2);
                       params(4,k) = mean(x)+colsStart-1;
                       params(5,k) = mean(y)+rowsStart-1;
                       k=k+1;
                   end
            end           
            
            if ~isempty(params)
                beads = magnificationCalibrationClass.testForNbeads(params,zDepth);
            end

        end
        
        %function to test whether image protion contains beads.
        function beads = testForNbeads(params,zDepth)
            x = params(4,:);
            y = params(5,:);
            k=1;
            beads = beadClass([],[],[],[]);
            while and(~isempty(x),~isempty(y));
                difX = abs(repmat(x(1),1,length(x))-x);
                difY = abs(repmat(y(1),1,length(y))-y);
                A = params(:,and(difX<5,difY<5));
                %This test the contours to see if there are more than 3
                %contours centred at a location, that there are at least 2
                %contours that have a radius greater or equal to 2. Additionally
                %need a mean of the stdev of contour radii to be <1,
                %testing for circular objects
                if size(A,2)>1 
                    beads(k) = beadClass(mean(A(3,:)),mean(A(4,:)),mean(A(5,:)),zDepth);
                    k=k+1;
                end
                params = params(:,~and(difX<5,difY<5));
                x = x(~and(difX<5,difY<5));
                y = y(~and(difX<5,difY<5));
            end
        end
        
        % function to get contour lines from matrix C. Taken from internet.
        function s = getcontourlines(c)
            sz = size(c,2);
            ii=1;
            jj=1;

            while ii<sz
                n = c(2,ii);
                s(jj).v = c(1,ii);
                s(jj).x = c(1,ii+1:ii+n);
                s(jj).y =  c(2,ii+1:ii+n);
                ii = ii +1 +n;
                jj=jj+1;
            end
        end
        
        %scatter points of beads on image
        function scatterPoints(beads,image)
            x = beads(1).centreX;
            y = beads(1).centreY;
            if size(beads,2)>1
                for i=2:size(beads,2)
                    x = horzcat(x,beads(i).centreX);
                    y = horzcat(y,beads(i).centreY);
                end
            end
            imagesc(image);
            hold on
            scatter(x,y,'r');
            hold off
            drawnow;
        end
            
            
    end
    
    %% private methods
    methods (Access = private)
        %test path type input and output path
        function path = testpath(obj,type)
            if strcmp(type,'fixed')==1
                    path = obj.fixedpath;
                elseif strcmp(type,'dynamic')==1
                    path = obj.dynamicpath;
                else
                    error('Please enter either fixed or dynamic');
            end
        end

    end               
    
end

