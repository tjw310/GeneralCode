classdef magnificationCalibrationClass < handle
    %Class and methods for the magnification calibration
    properties
        fixedpath
        dynamicpath
        imSt = []
        imFn = []
        DyMxIm = []
        StMxIm = []
        opticalCentre = []
        dmdy = []
        ERRdmdy = []
        dmtdy = 0.222;
        pxSz = 6.5e-3/20.33;
        
    end
    
    methods
        %% constructors and get/set
        %contructor
        function obj = magnificationCalibrationClass()
            path = uigetdir('Z:\','Select Folder of Images');
            obj.fixedpath = fullfile(path,'145mA');
            obj.dynamicpath = fullfile(path,'0-290mA');
        end
        
        %change ETL limits
        function ETLlimitchange(obj,low,high)
            obj.dynamicpath = fullfile(fileparts(obj.fixedpath),strcat(num2str(low),'-',num2str(high),'mA'));         
        end
        
        %change path constructor
        function changePath(obj)
            path = uigetdir('Z:\','Select Folder of Images');
            obj.fixedpath = fullfile(path,'145mA');
            obj.dynamicpath = fullfile(path,'0-290');
        end
        
        %% main functions
        % if it doesn't work or takes too long check the threshold value in
        % FastPeakFind.m
        function [fixed,mapped,magGrad,SR,beta,betaerr] = main(obj)
            fixed = obj.getFixedBeads;
            dynamic = obj.getDynamicBeads;
            mapped = obj.mapDepth(dynamic,fixed);
            obj.getOpticalCentre(mapped);
            magGrad = obj.scatterZmagGrad(mapped);
            SR = obj.scatterzSR(mapped);
            [beta,betaerr]=obj.beta(mapped);
            obj.dmdy = median(magGrad);
            obj.ERRdmdy = std(magGrad);
            obj.scatterAllPoints(mapped,'dynamic');
        end
        
        %Dynamic bead find
        function beadsOut = getDynamicBeads(obj,varargin)
            path = obj.testpath('dynamic');
            type = 'dynamic';
            if nargin>1
                path = varargin{1};
                type = path;
            end
            im_dir = dir(fullfile(path,'*.tiff'));
            if isempty(im_dir)
                im_dir = dir(fullfile(path,'*.tif'));
            end
            beadsOut = obj.trackBeads(type,1,size(im_dir,1));
           % beadsOut = obj.getFullTracks(beadsOut);
            %beadsOut = obj.subPixel(fullTracks,'dynamic');
        end
        
        %Static bead find
        function beadsOut = getFixedBeads(obj,varargin)
            % finds all fixed beads based on MIP of fixed stack, then finds
            % correct z-location looking for maximum of I-z profile at
            % given x,y location
            path = obj.testpath('fixed'); type = 'fixed';
            if nargin>1
                path = varargin{1}; %allows for different path
                type = path;
            end
            im_dir = dir(fullfile(path,'*.tiff'));
            if isempty(im_dir)
                im_dir = dir(fullfile(path,'*.tif'));
            end
            if isempty(obj.StMxIm)
                maxIm = loadimage(obj,type,1);
                size(im_dir)
                for i=2:size(im_dir,1)
                    maxIm = max(loadimage(obj,type,i),maxIm);
                    disp(i);
                end
                obj.StMxIm = maxIm;
            end
            
            [x,y] = magnificationCalibrationClass.findPeaks(obj.StMxIm,'show');
            
            mxbeads = [];
            thres = mean(obj.StMxIm(:))+2*std(obj.StMxIm(:));
            for i=1:length(x)
                if obj.StMxIm(y(i),x(i)) > thres
                    mxbeads = horzcat(mxbeads,beadClass(obj.StMxIm(y(i),x(i)),x(i),y(i),[]));
                end
            end
            
            beadsOut = [];
            for i=1:size(im_dir,1)
                [image,zDepth] = obj.loadimage(type,i);
                if i==1
                    for j=1:length(mxbeads)
                        beadsOut = horzcat(beadsOut,beadClass(image(mxbeads(j).centreY,mxbeads(j).centreX),mxbeads(j).centreX(1),mxbeads(j).centreY(1),zDepth));
                    end
                else
                    for j=1:length(beadsOut)
                        beadsOut(j).addPoint(beadClass(image(mxbeads(j).centreY,mxbeads(j).centreX),mxbeads(j).centreX(1),mxbeads(j).centreY(1),zDepth));
                    end
                end
                disp(i/size(im_dir,1)*100);
            end            
            beadsOut = obj.sortFixedBeads(beadsOut,type);          
        end
        
        %Mapping of dynamic traces to static positions
        function beadsOut = mapDepth(obj,dynBeads,fixedBeads)
            beadsOut = [];
            for i=1:length(dynBeads)
                for j=1:length(fixedBeads)
                    dx = abs(dynBeads(i).centreX-fixedBeads(j).centreX(1));
                    dy = abs(dynBeads(i).centreY-fixedBeads(j).centreY(1));
                    dr(j) = mean(sqrt(dx.^2+dy.^2));
                    mndr(j) = min(sqrt(dx.^2+dy.^2));
                end
                [~,lc] = min(dr);
                if mndr(lc) < 10
                    beadsOut = horzcat(beadsOut,mappedTrack(dynBeads(i),fixedBeads(lc).centreX(1),fixedBeads(lc).centreY(1),fixedBeads(lc).zDepth(1)));
                end
            end
        end
        
        %Find the optical centre of the mapped traces
        function getOpticalCentre(obj,mappedDynBeads)
            obj.scatterAllPoints(mappedDynBeads,'dynamic'); hold on;
            x = (1:size(obj.DyMxIm,2));
            k=1;
            for i=1:length(mappedDynBeads) 
                [a,b] = mappedDynBeads(i).fitxy;
                plot(x,a.*x+b); drawnow;
                if and(a~=inf,b~=inf) && and(~isnan(a),~isnan(b))
                    grad(k) = a;
                    inter(k)= b;
                    k=k+1;
                end
            end
            hold off;
               
            [x,y]=findClosestLineIntercept(grad,inter);
            obj.opticalCentre = [x,y];
            beadsOut = mappedDynBeads;
        end
        
        %Scatter relative magnification cone angle against static bead z
        %position and get magGrad
        function magGrad = scatterZmagGrad(obj,mappedDynBeads)
            for i=1:length(mappedDynBeads)
                magGrad(i) = mappedDynBeads(i).magGrad(obj.opticalCentre);
                z(i) = mappedDynBeads(i).statZ;
            end
            figure;
            scatter(z,magGrad,'x');
            xlabel('Static Z Position of Microspheres');
            ylabel('Magnification Gradient');
        end
        
        %Scatter scan range against static bead position. Additionally
        %outputs the relative maximum of the scan range, and half the scan
        %range (these should be approximately equal if everything is ok).
        function SR = scatterzSR(obj,mappedDynBeads)
            for i=1:length(mappedDynBeads)
                SR(i) = mappedDynBeads(i).getScanRange;
                z(i) = mappedDynBeads(i).statZ;
            end
            figure;
            scatter(z,SR,'x');
            xlabel('Static Z Position of Microspheres');
            ylabel('Scan Range');
        end
        
        %gets Beta, which is the source-detector distance.
        function [out,err] = beta(obj,mappedDynBeads)
            SR = obj.scatterzSR(mappedDynBeads);
            mz = obj.scatterZmagGrad(mappedDynBeads);
            mMZ = mean(mz); eMZ = std(mz);
            out = -1/mMZ;
            err = out*sqrt((eMZ/mMZ)^2);
        end
        
        function [xe,ze,xp,zp,y] = getMagandTrueBeadPosition(obj,images,fixed,angle)
            if angle==0
                image = images.projections(:,:,1);
                [x,z]=FastPeakFind(image.',1);
            elseif angle==180
                image = images.projections(:,:,201);
                [x,z]=FastPeakFind(image.',1);
            else
                error('Angle must be 0 or 180');       
            end
            
            figure; imagesc(image.'); hold on; scatter(x,z,'rx'); drawnow; hold off;

%             k=1;
%             
%             for i=1:length(fixed)
%                 r = sqrt((fixed(i).centreX-x).^2+(fixed(i).centreY-z).^2);
%                 [mn,lc] = min(r);
%                 fixed(i).realX = [];
%                 fixed(i).realY = [];
%                 if mn<50
%                     xe(k) = x(lc);
%                     ze(k) = z(lc);
%                     xp(k) = fixed(i).centreX;
%                     zp(k) = fixed(i).centreY;
%                     y(k) =  fixed(i).zDepth;
%                     locs(k) = i;
%                     k=k+1;
%                 end
%             end
            
            k=1;
            loopno = 1;
            for i=1:length(fixed)
                xf(i) = fixed(i).centreX;
                zf(i) = fixed(i).centreY;
                yf(i) = fixed(i).zDepth;
            end
            figure; scatter(x,z); hold on; scatter(xf,zf,'rx'); hold off; drawnow;
            while ~isempty(xf)
                r = sqrt((xf(1)-x).^2+(zf(1)-z).^2);
                [mn,lc] = min(r);
                if mn<50
                    xe(k) = x(lc);
                    ze(k) = z(lc);
                    xp(k) = xf(1);
                    zp(k) = zf(1);
                    y(k) =  yf(1);
                    locs(k) = loopno;
                    k=k+1;
                    
                    x = x((1:length(x))~=lc);
                    z = z((1:length(z))~=lc);
                end
                xf = xf(2:end);
                zf = zf(2:end);
                yf = yf(2:end);
                loopno = loopno+1;
            end

            figure; 
            %obj.scatterAllPoints(fixed,'fixed'); 
            hold on; scatter(xe,ze,[],locs); scatter(xp,zp,[],locs); 
            hold off;
        %figure;
            [e,~,X,Z,xe2,ze2] = findMagChanges(xp,xe,zp,ze,y,obj);
            
            images.e = e;

            for i=1:length(locs)
                idx = locs(i);
                fixed(idx).realX = X(i);
                fixed(idx).realY = Z(i);
            end

            dr = nanmedian(sqrt((xe2-xe).^2+(ze2-ze).^2));
            
            figure; imagesc(image.'); hold on; scatter(xe,ze,'ro');
            scatter(xp,zp,'gx'); scatter(xe2,ze2,'wx'); hold off; title(sprintf('dmdy = %.3f, medDist = %.3f',e,dr));
            end
            
        
        
        %% Sub-functions / Draw methods       
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
            if isempty(im_dir)
                im_dir = dir(fullfile(path,'*.tif'));
                imInfo = strsplit(im_dir(imageno).name,{'=','mA','.tif'});
            else
                im_dir(imageno).name;
                imInfo = strsplit(im_dir(imageno).name,{'=','mA','.tiff'});
            end
            zDepth =str2num(imInfo{2});
            image = double(imread(strcat(path,'\',char(cellstr(im_dir(imageno).name)))));
        end
        
        %load image from zDepth
        function image = loadimagefromZ(obj,type,zDepth)           
            path = obj.testpath(type);
            im_dir = dir(fullfile(path,'*.tiff'));
            if isempty(im_dir)
                im_dir = dir(fullfile(path,'*.tif'));
                for i=1:size(im_dir,1)
                    s = strsplit(im_dir(i).name,{'=','.tif'});
                    z(i) = str2num(s{2});
                end
            else
                for i=1:size(im_dir,1)
                    s = strsplit(im_dir(i).name,{'=','.tiff'});
                    z(i) = str2num(s{2});
                end
            end
            
            [~,lc]=min(abs(z-zDepth));
            z(lc);
            image = double(imread(strcat(path,'\',char(cellstr(im_dir(lc).name)))));
        end
        
        %get maximum intensity projection
        function mip(obj,type,varargin)
            path = obj.testpath(type);
            im_dir = dir(fullfile(path,'*.tiff'));
            if strcmp(type,'fixed')
                obj.StMxIm = obj.loadimage(type,1);
                for i=2:size(im_dir,1)
                    obj.StMxIm = max(obj.StMxIm,obj.loadimage(type,i));
                end
            elseif strcmp(type,'dynamic')
                obj.DyMxIm = obj.loadimage(type,1);
                for i=2:size(im_dir,1)
                    obj.DyMxIm = max(obj.DyMxIm,obj.loadimage(type,i));
                end
            end
            
            if ~isempty(varargin{1})
                subplot(2,1,1)
                imagesc(obj.StMxIm);
                subplot(2,1,2)
                imagesc(obj.DyMxIm);
            end               
        end
        
        %once have full bead traces, fit bead profiles to get sub-pixel
        %location.
        function beadsOut = subPixel(obj,beads,type)
            beadsOut = [];
            for i=1:length(beads)
                for j=1:length(beads(i).zDepth)
                    image = obj.loadimagefromZ(type,beads(i).zDepth(j));
                    [x,y,I]=obj.fitBead(image,beads(i).centreX(j),beads(i).centreY(j));
                    if j==1
                        beadsOut = horzcat(beadsOut,beadClass(I,x,y,beads(i).zDepth(j)));
                    else
                        beadsOut(i).addPoint(beadClass(I,x,y,beads(i).zDepth(j)));
                    end
                    disp(j/length(beads(i).zDepth)*100);
                end
                disp(i/length(beads)*100);
            end
        end
        
         %get beads from multiple images
        function tracked_Beads = trackBeads(obj,type,st,fn)
            obj.imSt = st;
            obj.imFn = fn;
            beads = [];
            if isempty(obj.DyMxIm)
                mxImflag =1;
            else
                mxImflag=0;
            end
            for i=st:fn
                if mxImflag==1
                    if i==st
                        obj.DyMxIm = obj.loadimage(type,st);
                    else
                        obj.DyMxIm = max(obj.DyMxIm,obj.loadimage(type,i));
                    end
                end
                
                nBeads = getBeads(obj,type,i);
                for j=1:size(nBeads,2)
                    if isempty(beads)
                        beads = nBeads(j);
                    else
                        for k=1:size(beads,2)
                            sz = size(sqrt((beads(k).centreX-nBeads(j).centreX).^2+(beads(k).centreY-nBeads(j).centreY).^2),2);
                            dif(k,1:sz) = sqrt((beads(k).centreX-nBeads(j).centreX).^2+(beads(k).centreY-nBeads(j).centreY).^2);
                        end
                        dif(dif==0)=NaN;
                        [r,c]=find(dif==min(dif(:)));
                        if length(r)>1
                            [r,lcmx] = max(r);
                            c = c(lcmx);
                        end
                        vl = dif(r,c);
                        
                        if vl<6
                            addPoint(beads(r),nBeads(j));
                        elseif vl>15
                            beads = horzcat(beads,nBeads(j));
                        %elseif vl<=40 && vl>=20 && dz>0.04
                         %   beads = horzcat(beads,nBeads(j));
                        end
                    end 
                end                  
                disp((i-st)/(fn-st)*100);
                obj.scatterAllPoints(beads,'dynamic');
                drawnow;
            end
            
            %Only include beads traces with more than 10 points.
            k=1;
            for j=1:size(beads,2)
                if length(beads(j).centreX)>10
                    tracked_Beads(k) = beads(j);
                    k=k+1;
                end
            end
            
            %display all points scattered on maximum intensity image
        end
        
        %get beads from image 
        function beads = getBeads(obj,type,imageno)
            [image,zDepth] = loadimage(obj,type,imageno);
            [x,y] = magnificationCalibrationClass.findPeaks(image);
            beads = [];         
            if ~isempty(x)
                for i=1:length(x)
                    beads = horzcat(beads,beadClass(image(y(i),x(i)),x(i),y(i),zDepth));
                end
            end
           % imagesc(image); hold on; scatter(x,y,'r'); hold off; drawnow;
        end
        
        %test bead tracks for full traces
        function fullTracks = getFullTracks(obj,tracked_Beads)
            k=1;
            while ~isempty(tracked_Beads)
                [flag,output]=tracked_Beads(1).testFullTrack;
                if flag==1
                    fullTracks(k) = output;
                    k=k+1;
                else
                    [flag,output,outIdx]=tracked_Beads(1).completeTrack(1,tracked_Beads);
                    if flag==1
                       tracked_Beads = tracked_Beads(1:length(tracked_Beads)~=outIdx);
                       fullTracks(k) = output;
                       k=k+1;
                    end
                end
                tracked_Beads = tracked_Beads(2:end);
            end   
        end
        
        %function to plot bead centre intensities vs zDepth
        function plotDepthIntensity(obj,beads)
            for i=1:size(beads,2)
                figure;
                scatter(beads(i).zDepth,beads(i).centIntensity);
            end
        end
        
        %function to scatter bead centres against maximum of image
        function scatterAllPoints(obj,allbeads,type,varargin)
            %figure;
            if ischar(type)
                if strcmp(type,'fixed')
                    imagesc(obj.StMxIm);
                else
                    imagesc(obj.DyMxIm);
                end
            else
                imagesc(type);
            end

            colormap(parula);
            hold on
            x=[];
            y=[];
            I=[];
            z=[];
            for i=1:size(allbeads,2)
                I = horzcat(I,allbeads(i).centIntensity);
                x = horzcat(x,allbeads(i).centreX);
                y = horzcat(y,allbeads(i).centreY);
                z = horzcat(z,allbeads(i).zDepth);
                if ~ischar(type)
                    x = x-(size(obj.DyMxIm,2)-size(type,2))/2;
                    y = y-(size(obj.DyMxIm,1)-size(type,1))/2;
                    X = allbeads(i).centreX-(size(obj.DyMxIm,2)-size(type,2))/2;
                    Y = allbeads(i).centreY-(size(obj.DyMxIm,1)-size(type,1))/2;
                else
                    X = allbeads(i).centreX;
                    Y = allbeads(i).centreY;
                end
                scatter(X,Y);
            end
            
            
            if ~isempty(obj.opticalCentre)
                scatter(obj.opticalCentre(1),obj.opticalCentre(2),'rx');
            end
            hold off
            drawnow;
            if nargin>3
                saveas(gcf(),fullfile(obj.dynamicpath,'tracks.fig'));
            end
%             figure;
%             imagesc(obj.mxIm);
%             hold on
%             scatter(x,y,[],z.*10000,'filled');
%             hold off
        end
        
        %scatter mapped tracks
        function scatterMapped(obj,mapTracks,varargin)
            imagesc(obj.DyMxIm);
            hold on
            x=[];
            y=[];
            I=[];
            z=[];
            X=[];
            Y=[];
            Z=[];
            for i=1:size(mapTracks,2)
                I = horzcat(I,mapTracks(i).centIntensity);
                x = horzcat(x,mapTracks(i).centreX);
                y = horzcat(y,mapTracks(i).centreY);
                z = horzcat(z,mapTracks(i).zDepth);
                X = horzcat(X,mapTracks(i).statX);
                Y = horzcat(Y,mapTracks(i).statY);
                Z = horzcat(Z,mapTracks(i).statZ);
            end
            scatter(x,y,'r');
            scatter(X,Y,'g');
            hold off
            drawnow;
            

%             figure;
%             hold on
%             scatter3(x,y,z,'r');
%             scatter3(X,Y,Z,'g');
%             hold off
        end
        
        %sort fixed beads function
        function beadsOut = sortFixedBeads(obj,beads,type)
            %This gets rid of any apparent beads that are not the correct
            %width. Correct width defined by total widths of all beads
            %found, so assumes the input beads contain a majority of true
            %beads. 
            beadsOut = [];
%            for i=1:length(beads)
%                I = beads(i).centIntensity;
%                [~,lc] = max(I);
%                locs = (lc-50:lc+50);
%                locs(locs<1)=1;
%                locs(locs>length(I))=length(I);
%                portion = I(locs);
%                [pks,~,w,~] = findpeaks(portion);
%                [~,idx]=sort(pks,'descend');
%                w = w(idx);
%                widths(i)=w(1);
%            end  

           for i=1:length(beads)
               %if widths(i)>mean(widths)-std(widths)
                   %[~,lc]=max(beads(i).centIntensity);
                   [~,lc,~,p] = findpeaks(beads(i).centIntensity);
                   [~,I]=sort(p,'descend');
                   lc = lc(I);
                   loc = lc(1);
                   drawnow;
                   wlocs = loc-30:loc+30;
                   wlocs(wlocs<1)=NaN;
                   wlocs(wlocs>length(beads(1).centIntensity)) = NaN;
                   wlocs = wlocs(~isnan(wlocs));
                   %wlocs(wlocs>length(beads(i).centIntensity)) = length(beads(i).centIntensity);
                   e = fit1DGaussian(beads(i).zDepth(wlocs),beads(i).centIntensity(wlocs));
                   if ~isnan(e(1))
                   [~,imlc] = min(abs(beads(i).zDepth-e(3)));
                   e(2)
                   %plot(beads(i).centIntensity); drawnow;
                   [image,zDepth]=obj.loadimage(type,imlc);
                   [x,y,I] = obj.fitBead(image,beads(i).centreX,beads(i).centreY);
                   [I,idx] = sort(I,'descend');
                   x=x(idx);
                   y=y(idx);
                   beadsOut = horzcat(beadsOut,beadClass(I(1),x(1),y(1),zDepth));
                   beadsOut(end).intensProfile = beads(i).centIntensity(wlocs);
                   beadsOut(end).zProfile = beads(i).zDepth(wlocs);
                   end
               %end
           end          
%            figure;
%            plot(widths);
%            hold on
%            plot(repmat(mean(widths)-std(widths),1,length(widths)));
%            hold off
        end
        
        %function to find index of bead nearest a position
        function idx = findBeadIdx(obj,allbeads,posX,posY)
            for i=1:size(allbeads,2)
                dr(i) = min((posX-allbeads(i).centreX).^2+(posY-allbeads(i).centreY).^2);
                [~,idx]=min(dr);
            end
        end
        
    end
    %% static methods
    methods (Static)
        %find points via FastPeakFind.m. This uses a local thresholding
        %method. Adapted from mathworks.
        function [x,y] = findPeaks(image,varargin)
            [x,y]=FastPeakFind(image,1);
            if nargin >1
                imagesc(image);
                hold on
                scatter(x,y,'r');
                hold off
                drawnow;
            end
        end
        
        %fits 2D gaussian to bead image to calculate sub-pixel central
        %location
        function [xOut,yOut,I,e] = fitBead(image,x,y)
            width = 15;
            coordsX = x-width:x+width;
            coordsX = coordsX(and(coordsX>=1,coordsX<=size(image,2)));
            coordsY = y-width:y+width;
            coordsY = coordsY(and(coordsY>=1,coordsY<=size(image,1)));
            
            [estimates,model,flag] = fit2DGaussian(image(coordsY,coordsX));
            
%            [~,fit]=model(estimates);           
%             subplot(1,2,1)
%             imagesc(image(coordsY,coordsX));
%             axis square
%             subplot(1,2,2)
%             imagesc(fit);
%             axis square
%             drawnow;
    
            e  = sqrt(1-(min(estimates(3:4))./max(estimates(3:4)))^2);  
            if flag==1
                xOut = estimates(5)+ coordsX(round(length(coordsX)/2));
                yOut = estimates(6)+ coordsY(round(length(coordsX)/2));
            else
                xOut = x;
                yOut = y;
            end
            I = estimates(1)+estimates(7);
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
                    path = type; %if not dynamic or fixed, assumes type is the path
            end
        end

    end               
    
end

