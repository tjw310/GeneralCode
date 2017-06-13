classdef beadClass < handle
    %Construct array for sequential images tracking central intensity and
    %location.
    
    properties
        centIntensity = []
        centreX = []
        centreY = []
        zDepth = []
        theta = []
        realX = []
        realY = []
        zProfile = []
        intensProfile = []
    end
    
    methods
        function obj = beadClass(varargin)
            if nargin==4
                obj.centIntensity = varargin{1};
                obj.centreX = varargin{2};
                obj.centreY = varargin{3};
                obj.zDepth = varargin{4};
            elseif nargin==0
                obj.centreX = [];
                obj.centreY = [];
                obj.centIntensity = [];
                obj.zDepth = [];
            elseif nargin==5
                obj.centIntensity = varargin{1};
                obj.centreX = varargin{2};
                obj.centreY = varargin{3};
                obj.zDepth = varargin{4};
                obj.theta = varargin{5};
            else
                error('Please enter valid properties');
            end
        end
        
        %add additional point
        function addPoint(obj,newObj)
            obj.centIntensity = horzcat(obj.centIntensity,newObj.centIntensity);
            obj.centreY = horzcat(obj.centreY,newObj.centreY);
            obj.centreX = horzcat(obj.centreX,newObj.centreX);
            obj.zDepth = horzcat(obj.zDepth,newObj.zDepth);
            if ~isempty(obj.theta)
                obj.theta = horzcat(obj.theta,newObj.theta);
            end
        end
        
        %merge tracks
        function objOut = mergeTracks(obj,newObj)
            objOut = beadClass();
            Inten = horzcat(obj.centIntensity,newObj.centIntensity);
            Y = horzcat(obj.centreY,newObj.centreY);
            X = horzcat(obj.centreX,newObj.centreX);
            Z = horzcat(obj.zDepth,newObj.zDepth);
            
            [objOut.zDepth,I] = sort(Z,'ascend');
            objOut.centreX = X(I);
            objOut.centreY = Y(I);
            objOut.centIntensity = Inten(I);
            
        end
        
        %plot bead centre intensities vs zDepth
        function plotZI(obj)
            plot(obj.zDepth,obj.centIntensity);
        end
        
        %scatter x,y location with filled intensity
        function scatterXYI(obj)
            scatter(obj.centreX,obj.centreY,[],obj.centIntensity,'filled');
        end

        %scatter x,y location with filled zDpeth
        function scatterXYZ(obj)
            scatter(obj.centreX,obj.centreY,[],obj.zDepth,'filled');
        end
        
        %3d scatter, x,y,z with I
        function scatterXYZI(obj)
            scatter3(obj.centreX,obj.centreY,obj.zDepth,[],obj.centIntensity,'filled');
        end
        
        %test for FULL bead tracks, minimum of 200 points, with two highly
        %prominent peaks
        function [exitFlag,varargout] = testFullTrack(obj)
            exitFlag =0;
            varargout{1} =[];
            if length(obj.zDepth)>40
                [~,trueLoc] = sort(obj.zDepth,'ascend');
                y = obj.centIntensity(trueLoc);
                %figure;
                findpeaks(y,'annotate','extents','MinPeakDistance',50);
                drawnow;
                [pks,locs,~,p]=findpeaks(y,'MinPeakDistance',25);
                [p,I] =sort(p,'descend');
                locs = locs(I);
                locs = trueLoc(locs);
                pks = pks(I);
                if p(1)>3*p(3) && p(2)>3*p(3) && pks(1)>3*min(obj.centIntensity) && pks(2)>3*min(obj.centIntensity) && (max(locs(1:2))-min(locs(1:2)))>180 && (max(locs(1:2))-min(locs(1:2)))<250
                    exitFlag = 1;
                    varargout{1} = splitTrack(obj,min(locs(1:2)),max(locs(1:2)));
                    %figure;
                    %plot(varargout{1}.centIntensity);
                end
            end
        end
        
        function objOut = splitTrack(obj,st,fn)
            objOut = beadClass(obj.centIntensity(st:fn),obj.centreX(st:fn),obj.centreY(st:fn),obj.zDepth(st:fn));
        end
        
        % if not a full track, attempts to find nearest other track, tests all tracks within 50 pixels 
        % of the end points, until has full track condition. If not returns
        % fail
        
        function [exitFlag,varargout] = completeTrack(obj,objIdx,allTracks)
            % varargout{1} - is the track itself
            % varargout{2} - is the additional track index.
            zSt = min(obj.zDepth);
            zEn = max(obj.zDepth);
            idx = find(obj.zDepth==min(obj.zDepth));
            spX = mean(obj.centreX(idx));
            spY = mean(obj.centreY(idx));
            idx = find(obj.zDepth==max(obj.zDepth));
            enX = mean(obj.centreX(idx));
            enY = mean(obj.centreY(idx));
            
            exitFlag = 0;
            varargout{1}=[];
            varargout{2} = [];
            k=1;
            ik=1;
            
            if size(allTracks,2)>1
                for i=1:size(allTracks,2)
                    if i==objIdx
                        rSt(i) = NaN;
                        rEn(i) = NaN;
                        dz(i) = NaN;
                    else      
                        rSt(i) = min(sqrt((spX-allTracks(i).centreX).^2+(spY-allTracks(i).centreY).^2));
                        rEn(i) = min(sqrt((enX-allTracks(i).centreX).^2+(enY-allTracks(i).centreY).^2));
                    end
                end
                
                [~,stNear] = find(rSt<100);
                [~,enNear] = find(rEn<100);
                lcs = [stNear,enNear];
                lcs = unique(lcs);
                
                goodflag=[];
                
                if ~isempty(lcs)
                    for i=1:length(lcs)
                        mergedTrack = mergeTracks(obj,allTracks(lcs(i)));
                        % if track less than 200 in length try and find
                        % another track, otherwise use 2 tracks
                        if length(mergedTrack.centreX)<200 && length(lcs)>1
                            m=1;
                            lcs2 = lcs(lcs~=lcs(i));
                            for j=1:length(lcs(lcs~=lcs(i)))
                                track = mergeTracks(mergedTrack,allTracks(lcs2(j)));
                                %figure;
                                %plot(track.centIntensity);
                                %drawnow;
                                [multiflag(j),clippedTrack] = testFullTrack(track);
                                if multiflag(j)==1
                                    tracks(m) = clippedTrack;
                                    idxes(m) = j;
                                    m=m+1;
                                end
                            end
                            if any(multiflag)
                                [clippedTrack,lc]=testForBestTrack(tracks);
                                goodTracks(k) = clippedTrack;
                                goodIdx(ik:ik+1)=lcs([i,idxes(lc)]);
                                k=k+1;
                                ik=ik+2;
                            end
                        else
                            [goodflag(i),clippedTrack] = testFullTrack(mergedTrack);
                            if goodflag(i)==1
                                goodTracks(k) = clippedTrack;
                                %plot(clippedTrack.centIntensity);
                                %drawnow;
                                goodIdx(ik) = lcs(i);
                                k=k+1;
                                ik=ik+1;
                            end
                        end
                    end
                end

                if any(goodflag)
                    exitFlag = 1;
                    if length(goodTracks)>=2
                        [varargout{1},lc] = testForBestTrack(goodTracks);
                        varargout{2} = goodIdx(lc);
                    else
                        varargout{1} = goodTracks;
                        varargout{2} = goodIdx;
                    end
                end 
            end
        end
        
        function [objOut,lc] = testForBestTrack(goodTracks)
            for i=1:length(goodTracks)
                p(i) = goodTracks(i).centIntensity(1)+goodTracks(i).centIntensity(2);
            end
            [~,lc] = max(p);
            objOut = goodTracks(lc);
        end
                
        function objOut = mergeTracksDirection(obj1,dir1,obj2,dir2)
            % dir = 0, leave array as is.
            % dir = 1, flip array back to front.
            objOut = beadClass();
            if dir1==1
                if dir2==1
                    objOut.centIntensity = horzcat(fliplr(obj1.centIntensity),fliplr(obj2.centIntensity));
                    objOut.centreY = horzcat(fliplr(obj1.centreY),fliplr(obj2.centreY));
                    objOut.centreX = horzcat(fliplr(obj1.centreX),fliplr(obj2.centreX));
                    objOut.zDepth = horzcat(fliplr(obj1.zDepth),fliplr(obj2.zDepth));
                else
                    objOut.centIntensity = horzcat(fliplr(obj1.centIntensity),obj2.centIntensity);
                    objOut.centreY = horzcat(fliplr(obj1.centreY),obj2.centreY);
                    objOut.centreX = horzcat(fliplr(obj1.centreX),obj2.centreX);
                    objOut.zDepth = horzcat(fliplr(obj1.zDepth),obj2.zDepth);
                end
            elseif dir2==1
                objOut.centIntensity = horzcat(obj1.centIntensity,fliplr(obj2.centIntensity));
                objOut.centreY = horzcat(obj1.centreY,fliplr(obj2.centreY));
                objOut.centreX = horzcat(obj1.centreX,fliplr(obj2.centreX));
                objOut.zDepth = horzcat(obj1.zDepth,fliplr(obj2.zDepth));
            else
                objOut.centIntensity = horzcat(obj1.centIntensity,obj2.centIntensity);
                objOut.centreY = horzcat(obj1.centreY,obj2.centreY);
                objOut.centreX = horzcat(obj1.centreX,obj2.centreX);
                objOut.zDepth = horzcat(obj1.zDepth,obj2.zDepth);
            end
        end
        




    end
    
end

