classdef mappedTrack < beadClass
    %Sub-class of beadClass determining corresponding static x,y,z position
    %of bead
    
    properties
        statX
        statY
        statZ
    end
    
    methods
        %constructor
        function obj = mappedTrack(beadObject,statX,statY,statZ)
            obj.statX = statX;
            obj.statY = statY;
            obj.statZ = statZ;
            obj.centreX = beadObject.centreX;
            obj.centreY = beadObject.centreY;
            obj.centIntensity = beadObject.centIntensity;
            obj.zDepth = beadObject.zDepth;
        end
        
        % fit linear relationship to x y
        function [grad,intercept,fit] = fitxy(obj)
            [estimates,~,~,fit]=fitLinear(obj.centreX,obj.centreY);
            grad = estimates(1);
            intercept = estimates(2);          
        end
        
        % get radius away from optical centre
        function out = r(obj,opticalCentre)
            out = sqrt((opticalCentre(1)-obj.centreX).^2+(opticalCentre(2)-obj.centreY).^2);
        end
        
        % get change in magnification
        function out = dm(obj,opticalCentre)
            rCent = sqrt((opticalCentre(1)-obj.statX).^2+(opticalCentre(2)-obj.statY).^2);
            out = obj.r(opticalCentre)./rCent;
        end
        
        % get change in z depth
        function out = dz(obj)
            out = -obj.statZ+obj.zDepth;
        end
        
        % get dm/dz
        function out = magGrad(obj,opticalCentre)
             [estimates,m] = fitLinear(obj.dz,obj.dm(opticalCentre));
             out = estimates(1);
             %[estimates,m] = fitPolynomial(obj.dz,obj.dm(opticalCentre));
             % out = estimates(6);
%             [~,f]=m(estimates);
%             scatter(obj.dz,obj.dm(opticalCentre));
%             hold on
%             plot(obj.dz,f);
%             hold off
        end

        % get scan range and maximum
        function [SR,mxSR] = getScanRange(obj)
            SR = mean(max(obj.dz))-mean(min(obj.dz));
        end
        % scatter change in magnification vs change in z depth
        function scatterzm(obj,opticalCentre)
            scatter(obj.dz,obj.dm(opticalCentre));
        end
        
        %3d scatter dx,dy,dz
        function scatterxyz(obj)
            scatter3(obj.dx,obj.dy,obj.dz,[],obj.centIntensity,'filled');
            xlabel('dx');
            ylabel('dy');
            zlabel('dz');
            colorbar();
            %clabel('Intensity');
        end
    end
    
end

