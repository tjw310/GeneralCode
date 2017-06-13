classdef coneRecon < handle
    %cone reconstruction class
    
    properties
        AoRcentreX %AoR centre in pixels in x
        AoRcentreY %AoR cnetre in pixels in y
        AoRmotion %axis of rotation motion as a function of angle in pixels
        AoRangle %angle of axis of rotation relative to the y axis
        AoRphiangle
        nProj
        nAngles
        nPx
        FOV
        dmtdy
        mtz0
        dmdy
        origOpticCentre %this is for the test object, as if overrites the optic centre if stage movement
        opticCentre
        interptype
        gpu
        filter
        projections
        coordsInit
        path
        filteredProj = []
        type
        rotOpticCent
        e
        stepsize
        stagemotion
        stmotX
        stmotZ
        stageangle
        piezomotion %total effective orbit (piezo position)
        fullpath
        DoF
        peaksX
        peaksZ
        peaksXsub
        peaksZsub        
    end
    
    methods
        %constructor
        function obj = coneRecon(varargin)
            obj.path = uigetdir();
            if nargin>0
                 obj.opticCentre = varargin{1};
                if nargin>1
                     obj.AoRcentreX = varargin{2};
                    if nargin>2
                         obj.AoRangle = varargin{3};
                    else
                        obj.AoRangle =0;
                    end
                else
                    obj.AoRcentreX = 0;
                    obj.AoRangle =0;
                end
            else
                obj.AoRcentreX = 0;
                obj.AoRangle =0;
            end
            obj.setParameters();
            obj.type = 'a'; %absorption type by default
        end
        %% trace predictions
        %use first and 2nd projections to predict parameters
        function predictParams(obj,fixed,trueX,trueZ)
            k = 100;
            t = 360-obj.theta;

            j=1;
            for i=1:length(fixed)              
                if ~isempty(fixed(i).realX)
                    cInit(1,j) = fixed(i).realX-size(obj.projections,2)/2;
                    cInit(3,j) = fixed(i).realY-size(obj.projections,1)/2;
                    zDepths(j) = fixed(i).zDepth;
                    beadno(j) = i;
                    j=j+1;
                end
            end
            ct=1;
            for i=1:size(cInit,2)
                cr = sqrt((cInit(1,i)-trueX(:,1)).^2+(cInit(3,i)-trueZ(:,1)).^2);
                [mini,loc]=min(cr);
                if mini<100
                    closeX(ct) = trueX(loc,1);
                    closeZ(ct) = trueZ(loc,1);
                    ct=ct+1;
                end
            end

           figure; imagesc(obj.projections(:,:,1).'); hold on; scatter(closeX+1024,closeZ+1024,'g');
           drawnow;

            %start_guesses = [obj.mtz0,obj.AoRcentre,obj.AoRcentre,0,-1,obj.opticCentre(1),obj.opticCentre(2)];
            start_guesses = [-50,10,-1,5];
            model = @fun;
            
            est = fminsearch(model,start_guesses)



            function [sse,x2,yR] = fun(params)
                mt = params(1);
                ARcentX = params(1);
                ARcentY = params(2);
                dy = obj.e(2);
                theta = params(3);
                gamma = params(4);
                
                op(1) = obj.opticCentre(1);
                op(2) = obj.opticCentre(2);

                %theta = 0;
                %gamma = 0;

                mt = obj.mtz0;
                %ARcentX = 0;
               % ARcentY = 0;
                u = cos(theta/180*pi)*sin(gamma/180*pi);
                v = sin(theta/180*pi)*sin(gamma/180*pi);
                w = cos(theta/180*pi);

                y = zDepths./obj.pxSz.*mt;

                x = (closeX-op(1)).*(1-dy.*y./mt*obj.pxSz)+op(1)-ARcentX;
                z = (closeZ-op(2)).*(1-dy.*y./mt*obj.pxSz)+op(2);
                
                y = y-ARcentY;

                xR = u*(u.*x+v.*y+w.*z)*(1-cos(t(k)/180*pi))+x*cos(t(k)/180*pi)+(-w.*y+v.*z)*sin(t(k)/180*pi);
                yR = v*(u.*x+v.*y+w.*z)*(1-cos(t(k)/180*pi))+y*cos(t(k)/180*pi)+(w.*x-u.*z)*sin(t(k)/180*pi);
                zR = w*(u.*x+v.*y+w.*z)*(1-cos(t(k)/180*pi))+z*cos(t(k)/180*pi)+(-v.*x+u.*y)*sin(t(k)/180*pi);

                %xR2 = x.*cos(t(k)/180*pi)-y.*sin(t(k)/180*pi);
                %yR2 = y.*cos(t(k)/180*pi)+x.*sin(t(k)/180*pi);

                xR = xR+ARcentX;
                yR = yR+ARcentY;

                x2 = (xR-op(1))./(1-dy.*yR./mt*obj.pxSz)+op(1);
                z2 = (zR-op(2))./(1-dy.*yR./mt*obj.pxSz)+op(2);
                
                c=1;
                for l=1:length(x2)
                    r = sqrt((x2(l)-trueX(:,k)).^2+(z2(l)-trueZ(:,k)).^2);
                    [mn,lc]=min(r);
                    if mn<100
                        s(c) = mn;
                        x3(c) = trueX(lc,k);
                        z3(c) = trueZ(lc,k);
                        c=c+1;
                    end
                end
                sse = sum(s.^4);
                %imagesc(obj.projections(:,:,k).'); 
                scatter(x2+1024,z2+1024,[],yR,'filled'); hold on; scatter(x3+1024,z3+1024,'rx'); hold off; drawnow;
            end
            
            %imagesc(obj.projections(:,:,k).'); hold on; scatter(x3+1024,z3+1024,'rx'); scatter(x2+1024,z2+1024,'g'); hold off; drawnow;
        end

        %fit single bead trace
        function fitSingleTracev2(obj,xrealArray,zrealArray,beadno,fixed)
            
            model = @funS;
            t = 360-obj.theta;
            
            xreal = xrealArray(beadno,:);
            zreal = zrealArray(beadno,:);
            
            yGuess = fixed(beadno).zDepth./obj.pxSz.*obj.mtz0;
            xGuess = (xreal(1)-obj.opticCentre(1)).*(1-obj.e(2).*yGuess./obj.mtz0*obj.pxSz)+obj.opticCentre(1);
            zGuess = (zreal(1)-obj.opticCentre(2)).*(1-obj.e(2).*yGuess./obj.mtz0*obj.pxSz)+obj.opticCentre(2);
            
            
            st = [obj.mtz0,0,0,xGuess,yGuess,zGuess,0,0];
            options = optimset('MaxFunEvals',2000*length(st));
            estimates = fminsearch(model,st,options)
            figure;
            funS(st);
            
            function sse = funS(params)            
                mt = params(1);
                %ARcentX = params(2);
                %ARcentY = params(3);

                dy = obj.e(2);
                phi = params(2);
                gamma = params(3);

                x = params(4);
                y = params(5);
                z = params(6);
                
                op(1) = obj.opticCentre(1);
                op(2) = obj.opticCentre(2);

                u = cos(phi/180*pi)*sin(gamma/180*pi);
                v = sin(phi/180*pi)*sin(gamma/180*pi);
                w = cos(phi/180*pi);
                
                for k=1:length(t)
                    model2 = @fun2;
                    st2 = [0,0];
                    ARestimates(:,k) = fminsearch(model2,st2,options)  ;       
                end
                
               
               scatter(xreal,zreal,'rx'); hold on; 
              % plot(x2,z2);
               scatter(x2,z2,'b');
               hold off; 
               drawnow;
               %pause(1);
               df = sqrt(20*(z2-zreal).^2+(x2-xreal).^2);
               %[df2,I] = sort(df,'descend');
               %df2 = df2(5:15);
               
               sse = nansum(df)
               
                function sse = fun2(params2)
                    ARcentX = params2(1);                    
                    ARcentY = params2(2);
                    yG = y-ARcentY;
                    xG = x-ARcentX;
                    xR = u*(u.*xG+v.*yG+w.*z)*(1-cos(t(k)/180*pi))+xG*cos(t(k)/180*pi)+(-w.*yG+v.*z)*sin(t(k)/180*pi)+ARcentX;
                    yR = v*(u.*xG+v.*yG+w.*z)*(1-cos(t(k)/180*pi))+yG*cos(t(k)/180*pi)+(w.*xG-u.*z)*sin(t(k)/180*pi)+ARcentY;
                    zR = w*(u.*xG+v.*yG+w.*z)*(1-cos(t(k)/180*pi))+z*cos(t(k)/180*pi)+(-v.*xG+u.*yG)*sin(t(k)/180*pi);
               
                    x2(k) = (xR-op(1))./(1-dy.*yR./mt*obj.pxSz)+op(1);
                    z2(k) = (zR-op(2))./(1-dy.*yR./mt*obj.pxSz)+op(2);
                    
                    sse = sqrt((x2(k)-xreal(k)).^2+(z2(k)-zreal(k)).^2);
                    
                    %scatter(xreal,zreal,'r');hold on; scatter(x2,z2,'b'); hold off; pause(0.05);
                end
            end
        end
        
        %fit bead traces
        function [est,dx,dz] = fitSingleTrace(obj,xrealArray,zrealArray,beadno,fixed)
            
            %model = @funS;
            t = 360-obj.theta;

            xreal = xrealArray(beadno,:);
            zreal = zrealArray(beadno,:);

%             st = [0,0];
%            options = optimset('MaxFunEvals',2e4*length(st),'MaxIter',2e4*length(st),'Tolx',1e-10,'Tolfun',1e-10);
%             estimates = fminsearch(model,st,options);
%             sse = model(estimates);
%             
%            [~,xplot,zplot,dx,dz,yR] = model(estimates);
            
            
            %figure; subplot(2,1,1); plot(dx);subplot(2,1,2); plot(dz);drawnow;
            
            model = @fun;
            yGuess = fixed(beadno).zDepth./obj.pxSz.*obj.mtz0;
            xGuess = (xreal(1)-obj.opticCentre(1)).*(1-obj.e(2).*yGuess./obj.mtz0*obj.pxSz)+obj.opticCentre(1);
            zGuess = (zreal(1)-obj.opticCentre(2)).*(1-obj.e(2).*yGuess./obj.mtz0*obj.pxSz)+obj.opticCentre(2);
            
            st = [xGuess,yGuess,zGuess];
            
            est = fminsearch(model,st);
            [~,xplot,zplot,dx,dz] = model(est);
            
            figure; subplot(2,1,1); plot(xplot,zplot); hold on; scatter(xreal,zreal,[],1:400); colormap(jet);
            scatter(obj.opticCentre(1),obj.opticCentre(2)); 
            hold off; title(num2str(beadno));
            subplot(2,1,2); plot(dx-mean(dx)); hold on; plot(dz-mean(dz)); hold off; drawnow;
            
            function [sse,x2,z2,dx,dz,yR] = funS(params)            
                %mt = params(1);
                mt = obj.mtz0;

                dy = obj.e(2);
                %phi = params(1);
                %gamma = params(2);
                phi = 0;
                gamma = 0;
                ARcentX = params(1);                    
                ARcentY = params(2);

                y = fixed(beadno).zDepth./obj.pxSz.*mt;
                x = (xreal(1)-obj.opticCentre(1)).*(1-dy.*y./mt*obj.pxSz)+obj.opticCentre(1);
                z = (zreal(1)-obj.opticCentre(2)).*(1-dy.*y./mt*obj.pxSz)+obj.opticCentre(2);
                
                
                
                op(1) = obj.opticCentre(1);
                op(2) = obj.opticCentre(2);

                u = cos(phi/180*pi)*sin(gamma/180*pi);
                v = sin(phi/180*pi)*sin(gamma/180*pi);
                w = cos(phi/180*pi);
                
                
                yG = y-ARcentY;
                xG = x-ARcentX;
                xR = u*(u.*xG+v.*yG+w.*z)*(1-cos(t/180*pi))+xG*cos(t/180*pi)+(-w.*yG+v.*z)*sin(t/180*pi)+ARcentX;
                yR = v*(u.*xG+v.*yG+w.*z)*(1-cos(t/180*pi))+yG*cos(t/180*pi)+(w.*xG-u.*z)*sin(t/180*pi)+ARcentY;
                zR = w*(u.*xG+v.*yG+w.*z)*(1-cos(t/180*pi))+z*cos(t/180*pi)+(-v.*xG+u.*yG)*sin(t/180*pi);
               
                x2 = (xR-op(1))./(1-dy.*yR./mt*obj.pxSz)+op(1);
                z2 = (zR-op(2))./(1-dy.*yR./mt*obj.pxSz)+op(2);
                
                
                
               % plot(z2); drawnow;
                
                xreal(or(xreal<=-1018,xreal>=1018)) = NaN;
                zreal(or(zreal<=-1018,zreal>=1018)) = NaN;
                
                dx = xreal-x2;
                dz = zreal-z2;
                                   
               %scatter(xreal,zreal,'rx'); hold on; 
              % plot(x2,z2);
              %scatter(x2,z2,'b');
              %hold off; 
              % drawnow;
               %pause(1);
               %df = sqrt(20*(z2-zreal).^2+(x2-xreal).^2);
               
               df = sqrt((z2-zreal).^2+(x2-xreal).^2);
               
               n = df>std(df);
               if n>0.1*length(df)
                   sse = 1000;
               else
                   sse = nanmean(df);
               end
               %[df2,I] = sort(df,'descend');
               %df2 = df2(5:15);
               
               
            end
            
            function [sse,x2,z2,dx,dz,yR] = fun(params)
                x1 = params(1);
                y1 = params(2);
                z1 = params(3);
                
                mt = obj.mtz0;
                dy = obj.e(2);
                
                xR = x1*cos(t/180*pi)-y1*sin(t/180*pi);
                yR = y1*cos(t/180*pi)+x1*sin(t/180*pi);
                zR = z1;
                
                x2 = (xR-obj.opticCentre(1))./(1-dy.*yR./mt*obj.pxSz)+obj.opticCentre(1);
                z2 = (zR-obj.opticCentre(2))./(1-dy.*yR./mt*obj.pxSz)+obj.opticCentre(2);
                
                x2(or(xreal<=-1018,xreal>=1018)) = NaN;
                z2(or(zreal<=-1018,zreal>=1018)) = NaN;
                
                dx = xreal-x2;
                dz = zreal-z2;
                
                sse = nansum(sqrt((z2-zreal).^2+(x2-xreal).^2));
            end
                
                
                
        end
        
        %get difference from trace
        function [tracks,dx,dz,xreal,zreal] = getTraceDif(obj,fixed,trueX,trueZ)
            j=1;
            for i=1:length(fixed)              
                if ~isempty(fixed(i).realX)
                    cInit(1,j) = fixed(i).realX-size(obj.projections,2)/2;
                    cInit(3,j) = fixed(i).realY-size(obj.projections,1)/2;
                    zDepths(j) = fixed(i).zDepth;
                    cInit(2,j) = fixed(i).zDepth./obj.pxSz.*obj.mtz0;
                    beadno(j) = i;
                    j=j+1;
                end
            end

            start_guesses = [obj.mtz0,obj.AoRcentre,obj.AoRangle,obj.e(2)];
            model = @fun;
           % estimates = fminsearch(model,start_guesses)
           
           [~,xreal,zreal] = fun(start_guesses);

            %obj.mtz0 = estimates(1);
            %obj.AoRcentre = estimates(2); obj.AoRangle = estimates(3);
            %obj.e(2) = estimates(4);
            
            function [sse,xreal,zreal] = fun(params)
                mt = params(1);
                cInit(2,:) = zDepths./obj.pxSz.*mt;
                ARcent = params(2);
                ARangle = params(3);
                dy = params(4);
                for k=1:length(obj.theta);
                    t = 360-obj.theta;
                    coords = coneRecon.T(cInit,[ARcent,0,0]);
                        if ARangle~=0
                            wVector = (1+tan(ARangle/180*pi)^2)^-0.5;
                            uVector = tan(ARangle/180*pi)*wVector;
                            coords = coneRecon.rotateComp(coords,uVector,0,wVector,t(k));
                        else
                            coords = coneRecon.rotateSimp(coords,t(k));
                        end
                        coords = coneRecon.T(coords,[-ARcent,0,0]);

                    tracks(1,:,k) = (coords(1,:)-obj.opticCentre(1))./(1-dy.*coords(2,:)./mt*obj.pxSz)+obj.opticCentre(1);
                    tracks(2,:,k) = coords(2,:);
                    tracks(3,:,k) = (coords(3,:)-obj.opticCentre(2))./(1-dy.*coords(2,:)./mt*obj.pxSz)+obj.opticCentre(2);

                    for jp=1:size(tracks,2)
                        dr = sqrt((tracks(1,jp,k)-trueX(:,k)).^2+(tracks(3,jp,k)-trueZ(:,k)).^2);
                        [mn,lc] = min(dr);
                        if mn<50
                            dx(jp,k) = tracks(1,jp,k)-trueX(lc,k);
                            dz(jp,k) = tracks(3,jp,k)-trueZ(lc,k);
                            xreal(jp,k) = trueX(lc,k);
                            zreal(jp,k) = trueZ(lc,k);
                            dR(jp,k) = dr(lc);
                        else
                            dx(jp,k)=1000;
                            dz(jp,k)=1000;
                            dR(jp,k) = 1000;
                            xreal(jp,k) = NaN;
                            zreal(jp,k) = NaN;
                        end
                    end 
                end
                
                sse = nansum(dR(:));
                %sse = nanmean(dR(:));
                dx(dx==1000)=NaN;
                dz(dz==1000)=NaN;
                dR(dR==1000)=NaN;
                %plot(dx.'); drawnow;
            end
            %figure; subplot(2,2,1); plot(dx.');
               %subplot(2,2,2); plot(dz.'); subplot(2,2,3); plot(dR.');
            
            cfill = squeeze(sqrt(tracks(1,:,:).^2+tracks(3,:,:).^2));
            cfill = squeeze(abs(tracks(1,:,:)));
            %figure; hold on;
            for i=1:size(dx,1)
                figure; hold on;
                %cfill = squeeze(atan2(tracks(1,i,:),tracks(3,i,:)));
                %scatter3(dx(i,:),tracks(2,i,:),dz(i,:),[],cfill);
                plot3(squeeze(tracks(1,i,:)),squeeze(tracks(2,i,:)),squeeze(tracks(3,i,:))); hold on;
                scatter3(xreal(i,:),squeeze(tracks(2,i,:)),zreal(i,:),'rx'); hold off;
%                 if mean(tracks(3,i,:))>0
%                      cl='r';
%                 else
%                     cl = 'b';
%                 end
%                 plot3(dx(i,:),squeeze(tracks(2,i,:)),dz(i,:),cl);
            end
hold off;

        end

        %create traces for real bead object based on using the fixed bead
        %calibration data
        function [tracks,beadno,lcs] = getTraces(obj,fixed)
            j=1;
            for i=1:length(fixed)              
                if ~isempty(fixed(i).realX)
                    disp(i);
                    cInit(1,j) = fixed(i).realX-size(obj.projections,2)/2;
                    cInit(3,j) = fixed(i).realY-size(obj.projections,1)/2;
                    zDepths(j) = fixed(i).zDepth;
                    cInit(2,j) = fixed(i).zDepth./obj.pxSz.*obj.mtz0;
                    beadno(j) = i;
                    j=j+1;
                end
            end

            figure;
            for k=1:length(obj.theta);
           % for k=1
                disp(k/length(obj.theta)*100);
                %t = max(obj.theta)-obj.theta;
                t = 360-obj.theta;
                coords = coneRecon.T(cInit,[obj.AoRcentre,0,0]);
                    if obj.AoRangle~=0
                        wVector = (1+tan(obj.AoRangle/180*pi)^2)^-0.5;
                        uVector = tan(obj.AoRangle/180*pi)*wVector;
                        coords = coneRecon.rotateComp(coords,uVector,0,wVector,t(k));
                    else
                        coords = coneRecon.rotateSimp(coords,t(k));
                    end
                    coords = coneRecon.T(coords,[-obj.AoRcentre,0,0]);
                    %coords(1,:) = coords(1,:)+obj.AoRmotion(k).*cos(obj.AoRangle/180*pi);
                    %coords(3,:) = coords(3,:)+obj.AoRmotion(k).*sin(obj.AoRangle/180*pi);
                
                tracks(1,:,k) = (coords(1,:)-obj.opticCentre(1))./(1-obj.e(2).*coords(2,:)./obj.mtz0*obj.pxSz)+obj.opticCentre(1);
                tracks(2,:,k) = coords(2,:);
                tracks(3,:,k) = (coords(3,:)-obj.opticCentre(2))./(1-obj.e(2).*coords(2,:)./obj.mtz0*obj.pxSz)+obj.opticCentre(2);
            
                if k==1
                imagesc(((1:size(obj.projections,1))-size(obj.projections,1)/2),...
                ((1:size(obj.projections,1))-size(obj.projections,1)/2),obj.projections(:,:,k).'); hold on;
                scatter(tracks(1,:,k),tracks(3,:,k),'wx'); hold off; drawnow; pause(0.2);
                end
            end
            
            figure; mip = squeeze(max(obj.projections,[],3));
            imagesc(((1:size(mip,1))-size(mip,1)/2),...
                ((1:size(mip,1))-size(mip,1)/2),mip.');
            hold on; for i=1:size(tracks,2)
                plot(squeeze(tracks(1,i,:)),squeeze(tracks(3,i,:)));
            end
            

            for i=1:size(tracks,2)
                [~,lcmxX(i)] = max(tracks(1,i,:));
                [~,lcmnX(i)] = min(tracks(1,i,:));
                [~,lcmxZ(i)] = max(tracks(3,i,:));
                [~,lcmnZ(i)] = min(tracks(3,i,:));
                scatter(tracks(1,i,lcmnZ(i)),tracks(3,i,lcmnZ(i)));
                scatter(tracks(1,i,lcmxZ(i)),tracks(3,i,lcmxZ(i)));
                scatter(tracks(1,i,lcmnX(i)),tracks(3,i,lcmnX(i)));
                scatter(tracks(1,i,lcmxX(i)),tracks(3,i,lcmxX(i)));                
            end
            
            hold off
            lcs = unique([lcmxX,lcmnX,lcmxZ,lcmnZ]); %unique for quad points
            lcs = unique([lcmxX,lcmnX]); %unique for two horizontal points

            lcs = vertcat(lcmxX,lcmnX); %gives bead information on which angles
            
            %[ estimates,flag ] = findTmag(obj,cInit,zDepths,lcs,beadno,mip)
        end

        %optimises the transverse magnification based on predicting bead
        %traces
        function [ estimates,flag ] = findTmag(obj,cInit,zDepths,lcs,beadno,mip)
        start_guesses = [obj.mtz0];
        model = @fun;
        options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
        t = repmat(360-obj.theta,2,1);
        thetaLcs = t(lcs);
            for row=1:size(lcs,1)
                for col=1:size(lcs,2)
                    [xra,zra] = FastPeakFind(obj.projections(:,:,lcs(row,col)).',1);
                    xr(1:length(xra),row,col) = xra-size(obj.projections,2)/2;
                    zr(1:length(zra),row,col) = zra-size(obj.projections,1)/2;
                end
            end

        [estimates,~,flag] = fminsearch(model,start_guesses,options);
        [~,mins] = model(estimates);
        figure; imagesc(mins); colorbar; drawnow;

            function [sse,mins] = fun(params)
                mt = params(1)
                ARcent = obj.AoRcentre;
                ARangle = obj.AoRangle;

                cInit(2,:) = zDepths/obj.pxSz*mt;

                for h=1:size(lcs,1)
                for k=1:size(lcs,2)

                cIn = cInit(:,beadno(k));
                coords = coneRecon.T(cIn,[ARcent,0,0]);
                    if ARangle~=0
                        wVector = (1+tan(ARangle/180*pi)^2)^-0.5;
                        uVector = tan(ARangle/180*pi)*wVector;
                        coords = coneRecon.rotateComp(coords,uVector,0,wVector,thetaLcs(h,k));
                    else
                        coords = coneRecon.rotateSimp(coords,thetaLcs(h,k));
                    end
                coords = coneRecon.T(coords,[-ARcent,0,0]);

                X = (coords(1,:)-obj.opticCentre(1))./(1-obj.e(2).*coords(2,:)./mt*obj.pxSz)+obj.opticCentre(1);
                Z = (coords(3,:)-obj.opticCentre(2))./(1-obj.e(2).*coords(2,:)./mt*obj.pxSz)+obj.opticCentre(2);


                imagesc(((1:size(obj.projections,1))-size(obj.projections,1)/2),...
                ((1:size(obj.projections,1))-size(obj.projections,1)/2),mip.'); hold on; scatter(xr(:,h,k),zr(:,h,k),'rx');
                scatter(X,Z,'g'); hold off; drawnow; pause(1);


                r = sqrt((xr(:,h,k)-X).^2+(zr(:,h,k)-Z).^2);
                mn=min(r);
                if mn<100
                    mins(h,k) = mn;
                else
                    mins(h,k) = NaN;
                end

                end
                end

            sse = nanmean(mins(:));
            end

        end
        
        %%
        %create test object and fill projections, tracks and coordsInit
        function tracks = testObject(obj,nBeads,rBeads,varargin)
            %figure;
            obj.opticCentre = obj.origOpticCentre;
            obj.filteredProj = [];
            obj.coordsInit = [];
            if nargin>3
            bool = 1;
            obj.coordsInit = varargin{1};
            else
            bool=0;
            end
            obj.projections = [];

            %x,y,z coords in object space, real size.
            sc = single((-obj.nPx/2:obj.nPx/2-1)*obj.pxSz/obj.mtz0).*obj.SubSampFc;
            dpx = obj.pxSz/obj.mtz0*obj.SubSampFc;
            [x,z] = meshgrid(gpuArray(single(sc)));
            y = gpuArray(single(sc));
            nslices = ceil(2*rBeads/dpx);

            tracks =[];
            %create object space and fill it will beads at coords
            XZ = gpuArray(single(zeros(obj.nPx,obj.nPx)));
            XY = gpuArray(single(zeros(obj.nPx,obj.nPx)));
            YZ = gpuArray(single(zeros(obj.nPx,obj.nPx)));
            plane = gpuArray(single(zeros(obj.nPx,obj.nPx)));
            XZbeads = gpuArray(single(zeros(obj.nPx,obj.nPx)));
            XYbeads = gpuArray(single(zeros(obj.nPx,obj.nPx)));
            YZbeads = gpuArray(single(zeros(obj.nPx,obj.nPx)));
           for k=1:length(obj.theta);
               % for k=1
                disp(k/length(obj.theta)*100);
                t = obj.theta;
                if k==1
                    if bool~=1
                        if nBeads>0
                            obj.coordsInit = coneRecon.getInitialObject(nBeads,obj); %get object initial coords
                        end
                    %obj.coordsInit = obj.addBeadAtOptCent(obj.coordsInit,obj.w);
                    end
                    %nBeads = nBeads+1;
                end
                    coords = coneRecon.T(obj.coordsInit,[obj.rAoRcentre,0,0]);
                    if obj.AoRangle==0 && obj.AoRphiangle==0
                        coords = coneRecon.rotateSimp(coords,t(k));                        
                    else
                        [uVector,vVector,wVector] = getUVW(obj.AoRangle,obj.AoRphiangle);
                        %wVector = (1+tan(-obj.AoRangle/180*pi)^2)^-0.5;
                        %uVector = tan(-obj.AoRangle/180*pi)*wVector;
                        coords = coneRecon.rotateComp(coords,uVector,vVector,wVector,t(k)); 
                    end
                    
                    

                    ypos(k,:) = coords(2,:);
                    rorbit(k,:) = sqrt(coords(1,:).^2+coords(2,:).^2+coords(3,:).^2);
                    
                
                    coords = coneRecon.T(coords,[-obj.rAoRcentre,0,0]);
                    coords(1,:) = coords(1,:)+obj.AoRmotion(k).*obj.pxSz/obj.mtz0/obj.SubSampFc;
                    coords(2,:) = coords(2,:)-1.3325*(obj.piezomotion(k)-mean(obj.piezomotion));
                    
                    %records tracks
                    [tracks(1,:,k),tracks(3,:,k)] = coneRecon.scaleCoords(coords(1,:)-obj.stagemotion(k),coords(2,:),coords(3,:),obj);
                    tracks(2,:,k) = coords(2,:);
                    tracks(1,:,k) = tracks(1,:,k) + obj.stagemotion(k);
                    
                    coords(1,:) = coords(1,:) - obj.stagemotion(k);
                    
                    for i=1:nBeads
                    [xz,xy,yz] = obj.getBeadProjection(coords(1,i),coords(2,i),coords(3,i),rBeads);
                    XZ = XZ + xz; XY = XY + xy; YZ = YZ + yz;
                    end
                    
                    
                    obj.opticCentre(1) = obj.opticCentre(1)+obj.stagemotion(k)/obj.pxSz*obj.mtz0/obj.SubSampFc.*cos(obj.stageangle/180*pi);
                    obj.opticCentre(2) = obj.opticCentre(2)+obj.stagemotion(k)/obj.pxSz*obj.mtz0/obj.SubSampFc.*sin(obj.stageangle/180*pi);
                   % imagesc(XZ); drawnow;
                    
                   %coords(3,:) = coords(3,:)+obj.AoRmotion(k).*obj.pxSz/obj.mtz0.*sin(obj.AoRangle/180*pi);
    
%                 for i=1:nBeads
%                     yidx = coords(2,i);
%                     [~,lc]=min(abs(y-yidx));
%                     N=0;
%                     [xCent,yCent,zCent]=coneRecon.coordsMagChange(coords(1,i),yidx,coords(3,i),obj);
% %                      if (xCent-rBeads)>(obj.pxSz/obj.mtz0*-obj.nPx/2) && (xCent+rBeads)<(obj.pxSz/obj.mtz0*(obj.nPx/2-1)) &&...
% %                              (yCent-rBeads)>(obj.pxSz/obj.mtz0*-obj.nPx/2) && (yCent+rBeads)<(obj.pxSz/obj.mtz0*(obj.nPx/2-1)) &&...
% %                                  (zCent-rBeads)>(obj.pxSz/obj.mtz0*-obj.nPx/2) && (zCent+rBeads)<(obj.pxSz/obj.mtz0*(obj.nPx/2-1))
%                         for j=-ceil(nslices/2):1:ceil(nslices/2) 
%                             n=0;
%                             if and(j+lc>0,(j+lc)<obj.nPx)
%                                 yC = y(j+lc);                           
%                                 [xidx,~,zidx]=coneRecon.coordsMagChange(xCent,-yC+y(lc),zCent,obj);
%                                 r2 = sqrt(rBeads^2-(dpx*abs(j))^2)./(1-obj.dmdy*yC);
%                                 r = sqrt((xidx-(x+obj.rStagemotionX(k))).^2+(zidx-(z+obj.rStagemotionZ(k))).^2);
%                                 bool1 = r<=r2;
%                                 if nnz(bool1)>0
%                                     plane(bool1)=plane(bool1)+1;
%                                     bool2 = and(r>r2,(r-r2)<(dpx));
%                                     if nnz(bool2)>0
%                                         plane(bool2) = plane(bool2)+(1-(r(bool2)-r2)./dpx);
%                                         n = sum(plane(:));
%                                     end
%                                     XZbeads = XZbeads+plane;
%                                 end
%                             N = N +n;
%                             XYbeads(j+lc,:) = sum(plane,1);
%                             YZbeads(:,j+lc) = sum(plane,2);
%                             plane(plane~=0)=0;
%                             end
%                         end
%                         %sum(XZbeads(:)./N)
%                         if strcmp(obj.type,'f')
%                          XY = XY + XYbeads./N;
%                          YZ = YZ + YZbeads./N;
%                          XZ = XZ + XZbeads./N;
%                         else
%                          XY = XY + XYbeads;
%                          YZ = YZ + YZbeads;
%                          XZ = XZ + XZbeads;
%                         end
%                      %end
%                      XYbeads(XYbeads~=0)=0;
%                      YZbeads(YZbeads~=0)=0;
%                      XZbeads(XZbeads~=0)=0;
%                 end
                
                
                

                %get xz projection for projection images
                obj.projections(:,:,k) = gather(XZ.');
                if rem(k,10)==0; obj.dispProjection(XZ,YZ,XY,tracks(:,:,k),k); end;
                
                XZ(XZ~=0)=0;
                XY(XY~=0)=0;
                YZ(YZ~=0)=0;
                
                obj.opticCentre = obj.origOpticCentre;
           end
           
%            [X,Z] = meshgrid(1:size(XZ,1));
%            for i=1:obj.nProj
%                im1 = obj.projections(:,:,i);
%                
%                X2 = (X-obj.opticCentre(1,1)/obj.SubSampFc+size(XZ,1)/2)./(1-obj.dmtdy.*obj.piezomotion(i))+obj.opticCentre(1,1)/obj.SubSampFc-size(XZ,1)/2;
%                Z2 = (Z-obj.opticCentre(1,2)/obj.SubSampFc+size(XZ,2)/2)./(1-obj.dmtdy.*obj.piezomotion(i))+obj.opticCentre(1,2)/obj.SubSampFc-size(XZ,2)/2;
%                obj.projections(:,:,i) = interp2(obj.projections(:,:,i).',X2,Z2).';
%                subplot(1,2,1); imagesc(X2-X);
%                subplot(1,2,2); imagesc(obj.projections(:,:,i)-im1); 
%                title(num2str(obj.piezomotion(i))); drawnow;
%            end
               
           %subplot(1,2,1); plot(ypos); hold on; plot(-obj.piezomotion.*obj.pxSz/obj.mtz0/obj.SubSampFc); hold off; subplot(1,2,2); plot(rorbit);

        end
        
        %reconstruct projections from FDK algorithm into volume given by
        %indexs (mnidx,mxidx)
        function recons = reconFV(obj,varargin)
            if nargin>1
                mnidx = varargin{1};
                mxidx = varargin{2};
            else
                mnidx=1;
                mxidx = obj.nPx;
            end
            if isempty(obj.filteredProj)
                obj.filteredProj = obj.filterProjections(obj.projections);
            end
            recons = obj.CTbackprojection(mnidx,mxidx);
            figure; subplot(2,2,1);
            imagesc((squeeze(max(recons,[],1))));
            subplot(2,2,2);
            imagesc((squeeze(max(recons,[],2))));
            subplot(2,2,3);
            imagesc((squeeze(max(recons,[],3))));
           
            figure;
            for i=1:size(recons,3)
                imagesc(recons(:,:,i));axis square;drawnow; pause(0.1);
            end

        end
        
        %reconstruct per slice and save to path
        function slice = reconPS(obj,mnidx,mxidx,varargin) 
            if isempty(obj.filteredProj)
                obj.filteredProj = obj.filterProjections(obj.projections);
            end
            %common parameters to send to reconstruction function
            [xx,yy] = meshgrid(obj.x,obj.y);
            t = obj.theta;
            if obj.gpu == 1
                xxS = gpuArray(xx-obj.AoRcentreX/obj.SubSampFc);
                yy = gpuArray(yy);
            else
                xxS = xx-obj.AoRcentreX/obj.SubSampFc;
            end
            
            for i=mnidx:mxidx
                tic
                slice = obj.CTbackprojectionSlicev2(i,xxS,yy,t);
                toc
                slice = gather(real(slice));
                slice = (slice+300).*10;
                %slice(isnan(slice))=0;
                %slice(slice<0)=0;
                if i==mnidx; figure; end; imagesc(slice);
                title(num2str(i)); axis square; colorbar;
                drawnow;
                slice = uint16(slice);
                if i==mnidx
                    if nargin>3
                        name = varargin{1};
      
                    else
                        name = 'reconSlices';
                    end
                    mkdir(obj.path,name);
                    p = fullfile(obj.path,name);
                end
                imwrite(slice,strcat(fullfile(p,num2str(i,'%05d')),'.tiff'));
            end
        end
        
        %reconstruct each central bead slice save to path
        function reconBeads(obj) 
            if isempty(obj.filteredProj)
                obj.filteredProj = obj.filterProjections(obj.projections);
            end
            figure;
            c = round((obj.coordsInit(3,:)./obj.pxSz*obj.mtz0))+obj.nPx/2;
            for i=1:length(c)
                idx = c(i);
                slice = obj.CTbackprojectionSlice(idx);
                slice = gather(uint16(slice));
                imagesc(slice);
                title(num2str(c(i)));
                drawnow;
                imwrite(slice,strcat(fullfile(obj.path,num2str(idx,'%05d')),'.tiff'));
            end
        end  
        
        %reconstruct from iradon
        function recons = reconiradon(obj,mnidx,mxidx)
            figure;
            for i=mnidx:mxidx
                sinogram = squeeze(obj.projections(:,i,:));
                sinogram = gpuArray(obj.shiftSinogram(sinogram));
                recons(:,:,i) = gather(iradon(sinogram,obj.theta,'ram-lak',1,size(sinogram,1)));
                disp(i);
                imagesc(recons(:,:,i));
            end
%             figure;
%             imagesc(squeeze(max(recons,[],2)));
%             figure;
%             imagesc(squeeze(sum(recons,2)));
        end
        
        %set parameters
        function setParameters(obj)
        %General objeters        
        obj.nProj = 400; %n.projections
        obj.nAngles = [0,360]; %angles covered
        obj.nPx = 512; %object space size in pixels;
        obj.FOV = 2160; %field of view in pixels

        % objeters from experimental data used to calcuate the cone height. 
        %object rotates around the z-axis. y-axis is depth
        %gradient of static magnification of fixed plane.
        obj.dmtdy = 0.2220;
        %transverse magnification at z=0
        obj.mtz0 = 20.3063;

        %relative magnification gradient of object
        obj.dmdy = -0.145;
        %obj.dmdy = -6;

        obj.interptype = 'linear'; % 'linear', 'nearest'
        obj.gpu = 1;
        obj.filter = 'ram-lak';
        
        obj.AoRmotion = zeros(1,obj.nProj);
        
        %clear filtered projections if filled
        obj.filteredProj = [];
        
        obj.stagemotion = zeros(1,400);
        obj.fullpath = zeros(1,400);
        obj.piezomotion = zeros(1,400);
        obj.stageangle = 0;
        
        obj.DoF = obj.FOV;
        
        end
         
        %get parameters from magnification calibration class
        function getParameters(obj,magClassObj)
            obj.origOpticCentre(1) = magClassObj.opticalCentre(1)-size(magClassObj.StMxIm,2)/2;
            obj.origOpticCentre(2) = magClassObj.opticalCentre(2)-size(magClassObj.StMxIm,1)/2;
            obj.dmdy = magClassObj.dmdy;
            obj.opticCentre = obj.origOpticCentre;
        end
        
        %load projection images from file (assume cropped to centre, if
        %cropped, varargin = nAngles either 180 or 360 (default is 360)
        function loadProjections(obj,varargin)
            obj.projections = [];
            obj.filteredProj = [];
            if nargin<=1   
                obj.path = uigetdir();
            end
            im_dir = dir(fullfile(obj.path,'*.tif'));
            if isempty(im_dir)
                im_dir = dir(fullfile(obj.path,'*.tiff'));
            end
            obj.nProj = length(im_dir);

            for i=1:obj.nProj
                if i==1
                    im_sz = size(single(imread(strcat(obj.path,'\',char(cellstr(im_dir(i).name))))).');
                    obj.projections = zeros(min(im_sz),min(im_sz),obj.nProj);
                end
                [im_name,matches] = strsplit(char(cellstr(im_dir(i).name)),{'y=','x=','z=','.tiff'});
                if size(im_name,2)>2
                    for k=1:length(matches)
                        switch matches{k}
                            case 'y='
                                obj.piezomotion(i) = str2double(im_name{k+1});
                            case 'x='
                                stageX(i) = str2double(im_name{k+1});
                            case 'z='
                                stageZ(i) = str2double(im_name{k+1});
                        end
                    end
                else
                    obj.piezomotion(i) = 0; stageX(i)=0; stageZ(i)=0; 
                end
                disp(i/obj.nProj*100);
                obj.nPx = min(im_sz);
                df = im_sz-obj.nPx;
                im = single(imread(strcat(obj.path,'\',char(cellstr(im_dir(i).name))))).';
                A = im(df(1)/2+1:end-df(1)/2,df(2)/2+1:end-df(2)/2);
                obj.projections(:,:,i) = A;
            end
            obj.nPx = size(obj.projections,1);
            obj.FOV = obj.nPx;
            obj.AoRcentreX = 0;
            obj.AoRangle = 0;
            
            
            stageX(isnan(stageX))=0; stageZ(isnan(stageZ))=0;
            stpX = abs(stageX-circshift(stageX,[0,1]));
            stpZ = abs(stageZ-circshift(stageZ,[0,1]));
            xStp = mean(stpX(stpX~=0)); xStp(isnan(xStp))=0;
            zStp = mean(stpZ(stpZ~=0)); zStp(isnan(zStp))=0;
            obj.stepsize = sqrt(xStp^2+zStp^2);
            if xStp~=0
                [~,f,~]=fitStageMotion(stageX,xStp);
                xFull = f-mean(f);
                xNorm = stageX-mean(f);
            else
                xFull = zeros(1,obj.nProj);
                xNorm = zeros(1,obj.nProj);
            end
            if zStp~=0
                [~,f,~]=fitStageMotion(stageZ,zStp);
                zFull = f-mean(f);
                zNorm = stageZ-mean(f);
            else
                zFull = zeros(1,obj.nProj);
                zNorm = zeros(1,obj.nProj);
            end
            
            obj.stmotX = xNorm; obj.stmotZ = zNorm;
%             obj.fullpath = xFull+zFull;
%             obj.fullpath = sqrt(xFull.^2+zFull.^2);
%             obj.fullpath(xFull<0) = -obj.fullpath(xFull<0);
%             obj.stagemotion = sqrt(xNorm.^2+zNorm.^2);
%             obj.stagemotion(xFull<0)  = -obj.stagemotion(xFull<0) ;
%             
             obj.stagemotion = round(xNorm/obj.stepsize).*obj.stepsize;
             
             s = obj.stagemotion; s(isnan(s))= 0; obj.stagemotion = s;
            
            if nargin>2
                obj.nAngles = [0,varargin{1}];
            else
                obj.nAngles = [0,360];
            end
            obj.projections = single(obj.projections);
            
            obj.stagemotion = obj.stagemotion(1:obj.nProj);
            obj.fullpath = obj.fullpath(1:obj.nProj);
            obj.piezomotion = obj.piezomotion(1:obj.nProj);
        end
        
        %reconstruct slice on optic centre. This should just be a fan beam
        %sinogram
        function reconOC(obj)
            sino = squeeze(obj.projections(:,round(size(obj.projections,2)/2+obj.opticCentre(2)),:));
            slice = obj.fanBeamrecon(sino);
        end
        
        function slice = fanBeamrecon(obj,sino,varargin) %optional inputs of mask
            if strcmp(obj.type,'f')
                fluor = 1; %fluorescence image
            elseif strcmp(obj.type,'a')
                fluor = 0; %absorption image
            else
                error('enter set type to f or a');
            end
            
            mask=0;
            if nargin>2
                if ~isempty(varargin)
                    mask = 1;
                end
            end
            u = obj.x-obj.opticCentre(1);
            for i=1:obj.nProj
                w = (obj.D^2+u.*sqrt(obj.R(i)^2-obj.D^2))./sqrt(obj.D^2+u.^2);
                sino(1:size(sino,1),i) = sino(:,i).*w';
            end

            filt_len = max(64,2^nextpow2(2*obj.nPx-1));
                    
            [ramp_kernel] = repmat(coneRecon.ramp_flat(obj.nPx),1,obj.nProj);
            
            lBefore = round((filt_len-obj.nPx)/2+obj.opticCentre(1));
            lAfter = filt_len-obj.nPx-lBefore;
            
            ramp_kernel= padarray(ramp_kernel,[lBefore,0],'pre');
            ramp_kernel = padarray(ramp_kernel,[lAfter,0],'post');
            
            sino = padarray(sino,[lBefore,0],'replicate','pre');
            sino = padarray(sino,[lAfter,0],'replicate','post');
            
            f = ifft(abs(fft(ramp_kernel)).*fft(sino));
            
            filterSino = f(lBefore+1:lBefore+obj.nPx,:);       
            
            t = obj.theta;
            [xx,yy] = meshgrid(single(obj.x),single(obj.y));
            v = single(zeros(obj.nPx,obj.nPx));
            if obj.gpu == 1
                yy = gpuArray(yy);     
                v2 = gpuArray(v-min(yy(:))+1);
                xxS = gpuArray(xx-obj.AoRcentreX);
            else
                xxS = xx-obj.AoRcentre;
            end
            
            slice = gpuArray(zeros(obj.nPx,obj.nPx,'single'));

            for i=1:obj.nProj
                if obj.dmdy>=0
                   beta = (360-t(i))/360*2*pi;
                else
                   beta = t(i)/360*2*pi;
                end

                u = obj.D.*(xxS*cos(beta)+yy*sin(beta)-obj.opticCentre(1)+obj.AoRposition(i))...
                   ./(xxS*sin(beta)-yy*cos(beta)+obj.D)-obj.AoRcentreX+obj.opticCentre(1);

                
                if obj.gpu == 1
                    u2 = gpuArray(u-xxS(1,1)+1);
                if fluor==1
                    if obj.dmdy>=0
                        dI = gpuArray((1+(xxS*sin(beta)-yy*cos(beta))*(-1*obj.dmdP)));
                    else
                        dI = gpuArray((1+(xxS*sin(beta)-yy*cos(beta))*obj.dmdP));
                    end
                end
                else
                    u2 = u-xxS(1,1)+1;
                    if fluor==1
                        if obj.dmdy>=0
                            dI = (1+(xxS*sin(beta)-yy*cos(beta))*(-1*obj.dmdP)).^2;
                        else
                            dI = (1+(xxS*sin(beta)-yy*cos(beta))*obj.dmdP).^2;
                        end
                    end
                end

                U = (xxS*sin(beta)-yy*cos(beta)+obj.D)./obj.D;
                proj = repmat(filterSino(:,i).',obj.nPx,1);
                %imagesc(proj); colorbar(); title(num2str(i)); drawnow;
                z = interp2(proj,u2,v2,obj.interptype);
                
                if fluor==1
                    recon = dI.*z./U.^2;
                else
                    recon = z./U.^2;
                end
                
                if mask==0
                    recon(isnan(recon))=0;
                end
 
                slice = slice + recon;
                
%                 if rem(i,5)==0
%                 subplot(2,2,1); imagesc(u2); colorbar(); axis square;
%                 title(num2str(i));
%                 subplot(2,2,2); imagesc(recon); colorbar(); axis square;
%                 subplot(2,2,3); imagesc(slice); colorbar(); axis square; drawnow;
%                 end
                                
            end
                slice(slice<0) = 0;
                slice = real(slice);
                figure;
                imagesc(slice); axis square; drawnow; 
                
        end    
       
        function projections = forwardProjectFanBeam(obj,slice)
            if strcmp(obj.type,'f')
                fluor = 1; %fluorescence image
            elseif strcmp(obj.type,'a')
                fluor = 0; %absorption image
            else
                error('enter type as f or a');
            end
            t = obj.theta;
            
            x = single(obj.x);
            y = single(obj.y);

            [xx,yy] = meshgrid(x,y);
            if obj.dmdy>=0
                f = single(repmat((1+obj.dmdP.*(-1*y)).',1,obj.nPx));
                A = single(repmat(((1+obj.dmdP.*(-1*y)).^2).',1,obj.nPx));
            else
                f = single(repmat((1+obj.dmdP.*(y)).',1,obj.nPx));
                A = single(repmat(((1+obj.dmdP.*(y)).^2).',1,obj.nPx));
            end
            if obj.gpu==1
                f = gpuArray(f);
                xx = gpuArray(xx);
                yy = gpuArray(yy);
            end
            for k=1:obj.nProj
                if obj.dmdy>=0
                    beta = (360-t(k))/360*2*pi;
                else
                    beta = t(k)/360*2*pi;
                end
                
                rotCent = obj.AoRposition(k);
                uu = ((xx-obj.opticCentre(1)).*f-rotCent+obj.opticCentre(1)).*cos(beta)-yy.*sin(beta)+obj.AoRcentre; 
                vv = ((xx-obj.opticCentre(1)).*f-rotCent+obj.opticCentre(1)).*sin(beta)+yy.*cos(beta);
                
                uu = uu-min(xx(:))+1;
                vv = vv-min(yy(:))+1;
                
                if fluor==1
                    r = A.*interp2(slice,uu,vv,obj.interptype);
                else
                    r = interp2(slice,uu,vv,obj.interptype);
                end
                projections(:,k) = nansum(r,1);
  
            end
        end
        
        function loop(obj)
            sino = squeeze(obj.projections(:,64,:));
            slice = obj.fanBeamrecon(sino,'mask');
            p = obj.forwardProjectFanBeam(slice);
            ob2 = obj.fanBeamrecon(p,'mask');
            p2 = obj.forwardProjectFanBeam(ob2);
            ob3 = obj.fanBeamrecon(p2,'mask');
            p3 = obj.forwardProjectFanBeam(ob3);
            ob4 = obj.fanBeamrecon(p3,'mask');
            p4 = obj.forwardProjectFanBeam(ob4);
            ob2 = ob2./max(ob2(:));
            ob3 = ob3./max(ob3(:));
            ob4 = ob4./max(ob4(:));

            figure; subplot(2,2,2); imagesc(ob3-ob2); axis square;
            subplot(2,2,1); imagesc(ob2-slice./max(slice(:))); axis square;
            subplot(2,2,3); imagesc(ob4-ob3); axis square;
            figure; subplot(2,2,1); imagesc(slice); axis square;
            subplot(2,2,2); imagesc(ob2); axis square;
            subplot(2,2,3); imagesc(ob3); axis square;
            subplot(2,2,4); imagesc(ob4); axis square;
            figure; subplot(2,2,3); imagesc(p3./max(p3(:))-p2./max(p2(:))); axis square;
            subplot(2,2,2); imagesc(p2./max(p2(:))-p./max(p(:))); axis square;
            subplot(2,2,4); imagesc(p4./max(p4(:))-p3./max(p3(:))); axis square;
            subplot(2,2,1); imagesc(p./max(p(:))-sino./max(sino(:))); axis square;
        end
        
        %find approx rotation angle and shift using MIP and correlation
        function getRotProj(obj)
            mip = squeeze(max(obj.projections,[],3));
            [obj.AoRangle,obj.AoRcentre] = coneRecon.findRotShift(mip);
            for i=1:size(obj.projections,3)
                rI = imrotate(padarray(obj.projections(:,:,i),[50,50],'replicate'),obj.AoRangle,'crop');
                obj.projections(:,:,i) = rI(51:end-50,51:end-50);
                disp(i);
            end
        end
        
        function findAoRcentre(obj,sino)
            %finds AoR centre by getting optic centre slice and fan beam
            %recon and forward project and compare sinograms. when have
            %difference between sinogram and forward project, ignore any
            %that are 3 sigma away         
            sLog = max(sino,[],2);
            idx = find(sLog>mean(sino(:)),1,'first');
            idx2 = find(sLog>mean(sino(:)),1,'last');
            st_guess = round((idx2+idx)/2-size(sino,1)/2);

            
            model = @fun;
            options = optimset('tolx',0.1,'tolfun',5e-4);
            offset = fminsearch(model,st_guess,options);
            
            function sse = fun(params)
                offset = params(1);
                obj.AoRcentre = offset;
                slice = obj.fanBeamrecon(sino,'mask');
                subplot(2,2,4); imagesc(slice); axis square; colorbar; drawnow;
                fp = gather(obj.forwardProjectFanBeam(slice));
                dif = abs(fp./max(fp(:))-sino./max(sino(:)));
                m = mean(dif(:));
                s = std(dif(:));
                dif(or(dif>m+3*s,dif<m-3*s)) = 0;
                sse = std(dif(:));              
            end
        end
            
        function out = getInitialObjectVolume(obj,rBeads)
            [x,y,z]=meshgrid(single(obj.rx));
            dpx = obj.pxSz/obj.mtz0;
            out = gpuArray(single(zeros(obj.nPx,obj.nPx,obj.nPx)));
            for i=1:size(obj.coordsInit,2)
                xC = obj.coordsInit(1,i);
                zC = obj.coordsInit(3,i);
                yC = obj.coordsInit(2,i);
                r = sqrt((xC-x).^2+(yC-y).^2+(zC-z).^2);
                if nnz(r<=rBeads)>0
                    out(r<=rBeads)=out(r<=rBeads)+1;
                    if nnz(and(r>rBeads,(r-rBeads)<dpx))>0
                        out(and(r>rBeads,(r-rBeads)<dpx)) = out(and(r>rBeads,(r-rBeads)<dpx))+(1-(r(and(r>rBeads,(r-rBeads)<dpx))-rBeads)./dpx);
                    end
                end
            end
        end

        function [X,Z,Xconc,Zconc] = findBeadsAllSlices(obj)
            Xconc = []; Zconc = [];
            %thres = nanmean(obj.projections(:))-sqrt(nanstd(obj.projections(:)));
            thres = nanmean(obj.projections(:))-nanstd(obj.projections(:))/5;
            for i=1:size(obj.projections,3)
                im = obj.projections(:,:,i);
                [z,x] = FastPeakFind(im,1);
                X(1:length(x),i) = x;
                Z(1:length(z),i) = z;
                Xconc = horzcat(Xconc,x);
                Zconc = horzcat(Zconc,z);
                disp(i);
                imagesc(im); hold on; scatter(z,x,'g'); hold off; drawnow;
            end
            X = X-size(obj.projections,2)/2;
            Z = Z-size(obj.projections,1)/2;
            Xconc = Xconc -size(obj.projections,2)/2;
            Zconc = Zconc -size(obj.projections,1)/2;
            
            obj.peaksX= X; obj.peaksZ = Z;
        end     
        
        function [XZ,XY,YZ] = getBeadProjection(obj,xc,yc,zc,radius)
            sigma = radius/sqrt(2*log(2));
            op = obj.RopticCentre;
            xo = op(1);
            zo = op(2);
            
            D = obj.rD;
            %% 2D test
%             [x,y] = meshgrid(obj.rx);      
%             x2 = (x-xo).*(1-y/D)+xo;
%             %I = exp(-((x-xc).^2+(y-yc).^2)./(2*sigma.^2)); 
%             %W = interp2(x,y,I,x2,y);   
%             weight = sqrt(x2.^2+D.^2); 
%             %W = W.*weight;
%             I2 = exp(-((x2-xc).^2+(y-yc).^2)./(2*sigma.^2));
%             I2 = I2.*weight;
%             XZ = zeros(obj.nPx); YZ = zeros(obj.nPx); XY = zeros(obj.nPx);
%               XZ(128,:) = nansum(I2,1);
%% 3D full array
%             [x,y,z] = meshgrid(obj.rx);
%             x2 = (x-xo).*(1-y/D)+xo;
%             z2 = (z-xo).*(1-y/D)+zo;
%             weight = sqrt(x2.^2+D.^2+z2.^2);
%             I2 = exp(-((x2-xc).^2+(y-yc).^2+(z2-zc).^2)./(2*sigma.^2));
%             I2 = I2.*weight;
%             XZ = squeeze(nansum(I2,1)).';
%             XY = squeeze(nansum(I2,3));
%             YZ = squeeze(nansum(I2,2)).';
%             XZ(isnan(XZ)) = 0; XY(isnan(XY)) = 0; YZ(isnan(YZ)) = 0;
%% 3D partial array
            dpx = obj.pxSz/obj.mtz0;
            nSlices = ceil(10*radius/dpx);
            [x,y,z] = meshgrid((-nSlices/2:nSlices/2-1).*dpx);
            if obj.dmdy<0;
                yc2= -yc; 
            else
                yc2 = yc; 
            end
            y2 = y+yc2;y = y+yc;
            xs = (xc-xo).*(1-yc2/D)+xo; zs = (zc-zo).*(1-yc2/D)+zo;
            x2 = (x-(xo-xc)).*(1-y2/D)+xo;
            z2 = (z-(zo-zc)).*(1-y2/D)+zo;
            weight = sqrt(x2.^2+D.^2+z2.^2);
            
            I2 = exp(-((x2-xs).^2+(y-yc).^2+(z2-zs).^2)./(2*sigma.^2));
            I2 = I2.*weight;
            
%             w = tukeywin(floor(obj.DoF/obj.SubSampFc),0.25);
%             w = padarray(w,[ceil((obj.FOV-obj.DoF)/(2*obj.SubSampFc)),0],'both');
%             w2 = interp1(obj.rx,w,squeeze(y(:,1,1)));
%             
%             g = repmat(w2,1,size(I2,2),size(I2,3));
%             I2 = I2.*g;
%             
%             subplot(2,2,3);
%             %imagesc(squeeze(y(:,1,1))+yc,squeeze(z(1,1,:))+zc,squeeze(nansum(I2,2)).');
%             imagesc(squeeze(nansum(I2,2)).');
%             xlabel('y /mm'); ylabel('z /mm'); subplot(2,2,2);
%             %imagesc(squeeze(x(1,:,1))+xc,squeeze(y(:,1,1))+yc,squeeze(nansum(I2,3)));
%             imagesc(squeeze(nansum(I2,3)));
%             xlabel('x /mm'); ylabel('y /mm'); subplot(2,2,1);
%             %imagesc(squeeze(x(1,:,1))+xc,squeeze(z(1,1,:))+zc,squeeze(nansum(I2,1)).');
%             imagesc(squeeze(nansum(I2,1)).');
%             xlabel('x /mm'); ylabel('z /mm'); %pause(0.4);
%             drawnow;
            
            xs2 = (xc-xo)./(1-yc2/D)+xo; zs2 = (zc-zo)./(1-yc2/D)+zo;
            
            [xr,zr] = meshgrid(squeeze(x(1,:,1))+xs2,squeeze(z(1,1,:))+zs2);
            [yr2,yr] = meshgrid(squeeze(y(:,1,1)),squeeze(y(:,1,1)));
            [xp,zp] = meshgrid(obj.rx); [yp2,yp] = meshgrid(obj.rx);
            XZ = interp2(xr,zr,squeeze(nansum(I2,1)).',xp,zp);
            XY = interp2(xr,yr,squeeze(nansum(I2,3)),xp,yp);
            YZ = interp2(yr2,zr,squeeze(nansum(I2,2)).',yp2,zp);
            XZ(isnan(XZ)) = 0; XY(isnan(XY)) = 0; YZ(isnan(YZ)) = 0;

            
        end
        
        function rotateProjections(obj,gamma)
            op = obj.opticCentre;
            obj.opticCentre(1) = op(1).*cos(gamma/180*pi)-op(2).*sin(gamma/180*pi);
            obj.opticCentre(2) = op(2).*cos(gamma/180*pi)+op(1).*sin(gamma/180*pi);
            for i=1:obj.nProj
                obj.projections(:,:,i) = imrotate(obj.projections(:,:,i),gamma,'crop');
            end
            obj.AoRangle = 0;
        end
        
        function [mtz0,motorMotion,AoRoffset] = findTmagandMotorMotion(obj)
            ymot = 1.33*(obj.piezomotion-mean(obj.piezomotion));
            xshift = circshift(ymot,[0,obj.nProj/4]);
            xdata = obj.peaksXsub*cos(obj.AoRangle/180*pi)-obj.peaksZsub*sin(obj.AoRangle/180*pi);
            zdata = obj.peaksZsub*cos(obj.AoRangle/180*pi)+obj.peaksXsub*sin(obj.AoRangle/180*pi);
            
            figure; scatter(xdata,zdata); drawnow;
           % xdata = obj.peaksXsub;

            start_guesses = [obj.mtz0,obj.AoRcentreX];
            model = @fun;
            model2 = @refine;
            options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));

            figure;
            estimates = fminsearch(model,start_guesses,options);
            estimates = fminsearch(model2,estimates,options)
            figure;
            [~,motorMotion] = model2(estimates); obj.AoRmotion = medfilt1(motorMotion,21,'omitnan','truncate');
      
            mtz0 = estimates(1); obj.mtz0 = mtz0;
            AoRoffset = estimates(2); obj.AoRcentreX = AoRoffset;
            %shiftpiezomotion = circshift(obj.piezomotion,[0,round(estimates(3))]);
           % obj.piezomotion = circshift(obj.piezomotion,[0,round(estimates(3))]);
            plot(motorMotion); hold on; plot(obj.AoRmotion); hold off;
            ylabel('Mean Motor Motion (pixels)'); xlabel('Proj No'); title(sprintf('With AoRoffset = %.2f',AoRoffset));
            savefig(gcf,strcat(obj.path,'/MotorMotion.fig'));

                function [sse,trueDf] = fun(params)
                        mt = params(1);
                        shift = params(2);
                        %dm = params(3);
                        dm=0;
                        %phaseshift = params(3);

                        %trueX = circshift(xshift./obj.pxSz.*mt,[0,round(phaseshift)])+shift;
                        
                        trueX = xshift./obj.pxSz.*mt+shift;
                        xfit = (trueX-obj.opticCentre(1))./(1-dm*ymot)+obj.opticCentre(1);
                        df = abs(xdata-xfit);
                        %df(df>50)=NaN;
                        sse = nansum(df.^2);

                        plot(xdata); hold on; plot(xfit); hold off; drawnow;
                        
                        dfScale = xdata-xfit;
                        trueDf = dfScale.*(1-dm*ymot);
                end
                
                function [sse,trueDf] = refine(params)
                        mt = params(1);
                        shift = params(2);
                       % dm = params(3);
                        dm=0;
                        %phaseshift = params(3);

                        %trueX = circshift(xshift./obj.pxSz.*mt,[0,round(phaseshift)])+shift;
                        trueX = xshift./obj.pxSz.*mt+shift;
                        xfit = (trueX-obj.opticCentre(1))./(1-dm*ymot)+obj.opticCentre(1);
                        df = abs(xdata-xfit);
                        df(df>50)=NaN;
                        sse = nansum(df.^2);

                        plot(xdata); hold on; plot(xfit); hold off; drawnow;
                        
                        dfScale = xdata-xfit;
                        trueDf = dfScale.*(1-dm*ymot);
                        trueDf(isnan(df)) = NaN;
                end
        end
        
        function findMagFromExtPoints(obj,amp,phs)
 
            xreal = obj.peaksXsub;
            model = @findMag;
            model2 = @refine;
            stG = [obj.mtz0,0,obj.AoRcentreX];
            estimates = fminsearch(model,stG);
            estimates = fminsearch(model2,estimates)
            [~,xpred] = model2(estimates);
            
            obj.AoRcentreX = estimates(3);
            
            figure; plot(xreal,'b'); hold on; plot(xpred,'r'); hold off;
            figure;plot(xreal-xpred); 
            
            function [sse,xpred] = findMag(params)
                a = params(1); %magnification
                b = params(2); %phase difference
                c = params(3); %axis offset
                
                xpred = a*1.333*amp*sin(obj.theta/180*pi+phs-2*pi/3+b-pi/2)/obj.pxSz+c;     
                
                %plot(xreal); hold on; plot(xpred); hold off; drawnow;
                %sse = sum(((xreal-xpred).*xreal).^2);
                
                [mx,lcmx] = max(xreal);
                [mn,lcmn] = min(xreal); 
                lcs = horzcat(lcmn-3:lcmn+3,lcmx-3:lcmx+3);
                lcs(lcs<1)=length(xreal)-lcs(lcs<1);
                lcs(lcs>length(xreal)) = lcs(lcs>length(xreal))-length(xreal);
                xrealpoints = xreal(lcs);
                xpredpoints = xpred(lcs);
                
                
               % plot(xrealpoints(1:floor(length(xrealpoints)/2))); hold on; plot(xpredpoints(1:floor(length(xrealpoints)/2))); hold off; drawnow;
               % plot(xrealpoints); hold on; plot(xpredpoints); hold off; drawnow;
                
                sse = sum((xrealpoints-xpredpoints).^2);
                
                
            end
            
            function [sse,xpred] = refine(params)
                a = params(1); %magnification
                b = params(2); %phase difference
                c = params(3); %axis offset
                
                xpred = a*1.333*amp*sin(obj.theta/180*pi+phs-2*pi/3+b-pi/2)/obj.pxSz+c;     
                
                %plot(xreal); hold on; plot(xpred); hold off; drawnow;
                %sse = sum(((xreal-xpred).*xreal).^2);
                
                [mx,lcmx] = max(xreal);
                [mn,lcmn] = min(xreal); 
                lcs = horzcat(lcmn-4:lcmn+4,lcmx-4:lcmx+4);
                lcs(lcs<1)=length(xreal)-lcs(lcs<1);
                lcs(lcs>length(xreal)) = lcs(lcs>length(xreal))-length(xreal);
                xrealpoints = xreal(lcs);
                xpredpoints = xpred(lcs);
                
                
               % plot(xrealpoints(1:floor(length(xrealpoints)/2))); hold on; plot(xpredpoints(1:floor(length(xrealpoints)/2))); hold off; drawnow;
               % plot(xrealpoints); hold on; plot(xpredpoints); hold off; drawnow;
                
               sse = nansum((xrealpoints-xpredpoints).^2);
               
               df = abs(xreal-xpred);
               df(df>50)=NaN;
               xpred(isnan(df)) = NaN;
                
                
            end
            
        end
                            
        function [Xout,Zout] = XZsubset(obj,xmn,xmx,zmn,zmx)
            %[x,z] = getTracesv2(obj);
            x = obj.peaksX; z = obj.peaksZ;
            %l = sum(and(and(x>xmn,x<xmx),and(z>zmn,z<zmx)),2);
           % [mx,lc] = max(l)
            for i=1:size(x,2)
                xa = x(:,i); za = z(:,i);
                xout(i,1:length(xa(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn))))) = xa(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn)));
                zout(i,1:length(za(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn))))) = za(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn)));
            end

            xout(xout==0)=NaN;
            zout(zout==0)=NaN;

            Xout = nanmean(xout,2).'; Zout = nanmean(zout,2).';
% 
%             X2 = Xconc; Z2 = Zconc;
%             X2(and(and(Z2<zmx,Z2>zmn),and(X2<xmx,X2>xmn))) = NaN;
%             Xout = Xconc(isnan(X2));
%             X2 = Xconc;
%             Z2(and(and(Z2<zmx,Z2>zmn),and(X2<xmx,X2>xmn))) = NaN;
%             Zout = Zconc(isnan(Z2));
%             figure; scatter(Xout,Zout,[],1:400,'filled'); 
%             disp(length(Xout));
            
            %Xout = x(lc,:); Zout = z(lc,:);
            obj.peaksXsub = Xout; obj.peaksZsub = Zout;

            
            e = fitLinear(Xout,Zout)
            theta = atan(e(1)); obj.AoRangle = theta/pi*180;

        end
            
        function rotProjections(obj)
            [X,Z]= meshgrid(obj.x,obj.z);           
            X2 = X*cos(obj.AoRangle/180*pi)-Z*sin(obj.AoRangle/180*pi)-X(1,1);
            Z2 = Z*cos(obj.AoRangle/180*pi)+X*sin(obj.AoRangle/180*pi)-Z(1,1);
            obj.opticCentre(1) = obj.opticCentre(1)*cos(obj.AoRangle/180*pi)-obj.opticCentre(2)*sin(obj.AoRangle/180*pi);
            obj.opticCentre(2) = obj.opticCentre(2)*cos(obj.AoRangle/180*pi)+obj.opticCentre(1)*sin(obj.AoRangle/180*pi);

            for i=1:obj.nProj
                disp(i/obj.nProj*100)
                single_rot_proj = single(interp2(obj.projections(:,:,i).',X2,Z2,obj.interptype));
                obj.projections(:,:,i) = single_rot_proj.';  
            end  
            obj.AoRangle = 0;
        end
        
        function scaleProjections(obj)
            [X,Z]= meshgrid(obj.x,obj.z);
            op(1) = obj.opticCentre(1);
            op(2) = obj.opticCentre(2);
            scFc = 1./(1-obj.dmtdy*(1.3330.*obj.piezomotion));
            for i=1:obj.nProj
                disp(i/obj.nProj*100)
                sc = scFc(i);
                X2 = (X-op(1)/obj.SubSampFc)*sc+op(1)/obj.SubSampFc-X(1,1);
                Z2 = (Z-op(2)/obj.SubSampFc)*sc+op(2)/obj.SubSampFc-Z(1,1);

                single_scaled_proj = single(interp2(obj.projections(:,:,i).',X2,Z2,obj.interptype));
                obj.projections(:,:,i) = single_scaled_proj.';  
            end
            obj.dmtdy = 0;
        end


    end
    
    methods %methods to get alternate properties and plot functions     
        %% methods to call necessary attributes
        function getProperties(obj,obj2)
            %get properties from another coneRecon object
            obj.AoRcentreX = obj2.AoRcentreX;
            obj.AoRcentreY = obj2.AoRcentreY;
            obj.AoRmotion=obj2.AoRmotion;
            obj.AoRangle =obj2.AoRangle ;
            obj.AoRphiangle= obj2.AoRphiangle;
            obj.nProj=obj2.nProj;
            obj.nAngles=obj2.nAngles;
            obj.FOV=obj2.FOV;
            obj.dmtdy=obj2.dmtdy;
            obj.mtz0=obj2.mtz0;
            obj.dmdy=obj2.dmdy;
            obj.origOpticCentre=obj2.origOpticCentre;
            obj.opticCentre=obj2.opticCentre;
            obj.interptype=obj2.interptype;
            obj.gpu=obj2.gpu;
            obj.stagemotion=obj2.stagemotion;
            obj.piezomotion=obj2.piezomotion;
            obj.DoF=obj2.DoF;
        end
        
        function out = rx(obj) %real x scale
            out = obj.x.*obj.pxSz./obj.mtz0;
        end
        
        function out = ry(obj)%real y scale
            out = obj.y.*obj.pxSz./obj.mtz0;
        end
        
        function out = rz(obj)%real z scale
            out = obj.z.*obj.pxSz./obj.mtz0;
        end
        
        function out = theta(obj) %calculate theta angles in degrees
            out = obj.nAngles(1):(obj.nAngles(2)-obj.nAngles(1))/obj.nProj:obj.nAngles(2)-(obj.nAngles(2)-obj.nAngles(1))/obj.nProj;
        end
        
        function out = SubSampFc(obj) %subsampling factor;
            out = obj.FOV/obj.nPx;
        end
        
        function out = w(obj)
            out = obj.nPx*obj.pxSz/obj.mtz0; %real size in mm of object space
        end
        
        function out = SR(obj) %scan range is equal to the real width of object
            out = obj.w;
        end
        
        function out = dmdP(obj)
            if obj.dmdy<0
            out = -(1/obj.rD).*obj.pxSz./obj.mtz0; %relative mag change per z pixel (z goes about 0)
            else
                out = (1/obj.rD).*obj.pxSz./obj.mtz0;
            end
        end
        
        function out = rD(obj) %distance from source to detector in mm
            out = abs(1/(-obj.dmtdy-obj.dmdy));
            out = abs(1/-obj.dmdy)+obj.AoRcentreY;
            out = abs(1/-obj.dmdy);
        end
        
        function out = D(obj) %distance from source to detector in pixels
            %out = abs(-1./obj.dmdP);
            % this includes the scaling factor for the magnification
            % change. It assumes that z/D is small.
            out = abs(obj.rD./obj.pxSz.*obj.mtz0);
        end
        
        function out = R(obj,projNo)   %distance from source to AOR in pixels
            L = obj.opticCentre(1).*cos(obj.AoRangle/180*pi)-obj.opticCentre(2).*sin(obj.AoRangle/180*pi);
            out = sqrt((L-obj.AoRposition(projNo))^2+obj.D^2); %pixels
        end
        
        function out = coneAngle(obj)
            out = atan2(sqrt((abs(obj.opticCentre(1))+obj.nPx/2)^2+(abs(obj.opticCentre(2))+obj.nPx/2)^2),obj.D)/pi*180;
        end
        
        function out = x(obj)
            out = (-obj.nPx/2:obj.nPx/2-1); % scale of object x in pixels
        end
        
        function out = y(obj)
            out = (-obj.nPx/2:obj.nPx/2-1); % scale of object y in pixels
        end
        
        function out = z(obj)
            out = (-obj.nPx/2:obj.nPx/2-1); % scale of object z in pixels
        end
        
        function out = u(obj) %detector system. x coordinates offset by distance. in pixels
            out = obj.x - obj.opticCentre(1);
        end
        
        function out = v(obj) %detector y coord system. (same as z). in pixels
            out = obj.z;
        end
        
        function out = RopticCentre(obj) %real distance optic centre in mm
            out = obj.opticCentre*6.5e-3./obj.mtz0;
        end
        
        function out = pxSz(obj) %pxSz in mm, if binned leads to an effective size
            out = 6.5e-3*obj.SubSampFc;
        end
        
        function out = rAoRcentre(obj)
            out = obj.AoRcentreX.*obj.pxSz/obj.mtz0/obj.SubSampFc;
        end %real AoR centre position in mm
        
        function out = rAoRposition(obj,angle)
            out = obj.AoRposition(angle).*obj.pxSz/obj.mtz0;
        end %gets current real AoR position in mm
        
        function out = stagemotionX(obj,angle)
            out = obj.stagemotion(angle).*cos(obj.stageangle/180*pi);
        end %gets current real AoR position in mm
        
        function out = stagemotionZ(obj,angle)
            out = obj.stagemotion(angle).*sin(obj.stageangle/180*pi);
        end %gets current real AoR position in mm

          
    %% plot functions
        function scatterCoords(obj,varargin) %scatter initial test object coordinates
            if nargin>1
                c = varargin{1};
                projNo = varargin{2};
            else
                c = obj.coordsInit;
                projNo = 1;
            end
            scatter3(c(1,:),c(2,:),c(3,:),'filled');
            xlim([min([obj.rx]),max([obj.rx])]);
            ylim([min([obj.ry]),max([obj.ry])]);
            zlim([min([obj.rz]),max([obj.rz])]);
            xlabel('x'); ylabel('y'); zlabel('z');
            hold on
            scatter3(repmat(obj.rAoRposition(projNo),1,obj.nPx),zeros(1,obj.nPx),obj.rz); hold off;
            drawnow;
        end
        
        function plotTracks(obj) %plots the tracks
            figure;
            for l=1:size(obj.tracks,2)   
                subplot(1,3,3)
                plot(squeeze(obj.tracks(1,l,:))/obj.w*obj.nPx,squeeze(obj.tracks(2,l,:))/obj.w*obj.nPx); hold on;
                if l==1
                    scatter(obj.opticCentre(1),0,'rx');
                    scatter(obj.AoRcentre,0);
                end
                xlabel('x');ylabel('y');
                subplot(1,3,2)
                plot(squeeze(obj.tracks(2,l,:))/obj.w*obj.nPx,squeeze(obj.tracks(3,l,:))/obj.w*obj.nPx); hold on;
                if l==1
                    scatter(0,obj.opticCentre(2),'rx');
                end
                xlabel('y'); ylabel('z');
                subplot(1,3,1)
                plot(squeeze(obj.tracks(1,l,:))/obj.w*obj.nPx,squeeze(obj.tracks(3,l,:))/obj.w*obj.nPx); hold on;
                if l==1
                    scatter(obj.opticCentre(1),obj.opticCentre(2),'rx');
                    
                    plot(repmat(obj.AoRcentre,1,obj.nPx),obj.z);
                end
                xlabel('x'); ylabel('z');
            end
            hold off;
        end
        
        function dispProjection(obj,XZ,YZ,XY,coords,projNo)
            op = obj.RopticCentre;
            subplot(2,2,1)
            imagesc(obj.rx+obj.stagemotionX(projNo),obj.rz+obj.stagemotionZ(projNo),XZ); xlabel('x /mm');ylabel('z /mm'); axis square;
            set(gca,'ydir','normal');
            hold on
            scatter(op(1),op(2),'go');
            y = tan(-obj.AoRangle/180*pi).*(obj.rx-obj.rAoRposition(projNo));
            plot(y,obj.rz);
            scatter(coords(1,:),coords(3,:),'rx');
            hold off
            subplot(2,2,2)
            imagesc(obj.rx+obj.stagemotionX(projNo),obj.ry,XY); axis square; set(gca,'ydir','normal'); xlabel('x /mm');ylabel('y /mm');
            hold on
            scatter(op(1),0,'go');
            scatter(obj.rAoRposition(projNo),0);
            scatter(coords(1,:),coords(2,:),'rx');
            hold off
            subplot(2,2,3)
            imagesc(obj.ry,obj.rz+obj.stagemotionZ(projNo),YZ); axis square; set(gca,'ydir','normal'); xlabel('y /mm');ylabel('z /mm');
            hold on
            scatter(0,op(2),'go');
            scatter(coords(2,:),coords(3,:),'rx');
            hold off
            drawnow;
        end
    end
    
    methods (Access = private)

        function out = AoRposition(obj,projNo) %finds AoR shift as a function of angle in degrees
            out = (obj.AoRcentreX + obj.AoRmotion(projNo));
        end
        
        function filt = generateFilter(obj)
            %FILTERING with x,y,z centred on object with u offset on optic axis.    
            filt_len = max(64,2^nextpow2(2*obj.nPx-1));
                    
            [ramp_kernel] = repmat(coneRecon.ramp_flat(obj.nPx),1,obj.nPx);

            lBefore = round((filt_len-obj.nPx)/2+obj.opticCentre(1));
            lAfter = filt_len-obj.nPx-lBefore;
            
            ramp_kernel= padarray(ramp_kernel,[lBefore,0],'pre');
            ramp_kernel = padarray(ramp_kernel,[lAfter,0],'post');
            
            filt = gpuArray(abs(fft(ramp_kernel)));
            
            filt = imrotate(filt,obj.AoRangle,'crop');
        end
        
        function [ proj ] = filterProjections(obj,proj)
            %FILTERING with x,y,z centred on object with u offset on optic axis. 
            u = linspace(-obj.nPx/2,obj.nPx/2,obj.nPx);
            v = linspace(-obj.nPx/2,obj.nPx/2,obj.nPx);
            [uuS,vvS] = meshgrid(u,v);
            
            if obj.dmdy<0
                gamma = -obj.AoRangle/180*pi;
            else
                gamma = -obj.AoRangle/180*pi;
            end
                        
            op = obj.opticCentre/obj.SubSampFc;
            
            filt_len = max(64,2^nextpow2(2*obj.nPx-1));
            
            mxfreq = 1/(obj.pxSz/obj.mtz0);
            mxNAfreq = 2*0.5/515e-6;
            fract = (mxNAfreq/mxfreq/2);
            %fract = 1;
            len = fract*filt_len-rem(fract*filt_len,2);
            
            window = tukeywin(len,0.01);
            window = repmat(padarray(window,[(filt_len-len)/2,0]),1,obj.nPx);

            
            [ramp_kernel] = repmat(coneRecon.ramp_flat(obj.nPx),1,obj.nPx);

            lBefore = round((filt_len-obj.nPx)/2+op(1));
            lAfter = filt_len-obj.nPx-lBefore;
            
            ramp_kernel= padarray(ramp_kernel,[lBefore,0],'pre');
            ramp_kernel = padarray(ramp_kernel,[lAfter,0],'post');
            
            %window = padarray(window,[lBefore,0],'pre');
            %window = padarray(window,[lAfter,0],'post');
            window = circshift(window,[filt_len/2,0]);
            
            filt = gpuArray(abs(fft(ramp_kernel)));
            
            figure;
            subplot(1,2,1)
            imagesc(filt);
            
            filt = filt.*window;
            subplot(1,2,2)
            imagesc(filt);
            drawnow;
            figure;
            
            [xxS,yyS] = meshgrid(obj.x);
            
            xx = xxS.*cos(gamma)-yyS.*sin(gamma);
            yy = yyS.*cos(gamma)+xxS.*sin(gamma);
            
            xx = xx-xxS(1,1)+1;
            yy = yy-yyS(1,1)+1;
            
            xx2 = xxS.*cos(gamma)+yyS.*sin(gamma);
            yy2 = yyS.*cos(gamma)-xxS.*sin(gamma);
            
            xx2 = xx2-xxS(1,1)+1;
            yy2 = yy2-yyS(1,1)+1;
            
            for i=1:obj.nProj
                disp(i/obj.nProj*100);
                s = obj.stagemotion(i);
                op(1) = op(1)*cos(gamma)-op(2)*sin(gamma);
                op(2) = op(2)*cos(gamma)+op(1)*sin(gamma);
                 subplot(1,3,1); imagesc(proj(:,:,i));
                
                p = proj(:,:,i);
                
                p = interp2(p,xx,yy);
                p(isnan(p))=0;
                
                uu = uuS-op(1);
                vv = vvS-op(2);
                w = (1-uu.*op(1)/obj.SubSampFc./(obj.D^2+vv.^2))./sqrt(obj.D^2+uu.^2 + vv.^2);
                p = p.*w';
   
                fproj = padarray(p,[lBefore,0],'replicate','pre');
                fproj = padarray(fproj,[lAfter,0],'replicate','post');
                fproj = gpuArray(fproj);

                fproj = fft(fproj);

                fproj = fproj.*filt;
                subplot(1,3,3); imagesc(filt);

                fproj = real(ifft(fproj));
  
               proj(:,:,i) = gather(fproj(lBefore+1:lBefore+obj.nPx,:));
               
               p2 = interp2(proj(:,:,i),xx2,yy2);
               p2(isnan(p2))=0;
               %p2 = imfilter(p2,fspecial('laplacian',.0001),'replicate');
               proj(:,:,i) = p2;
               
               subplot(1,3,2); imagesc(proj(:,:,i)); drawnow;
            end
        end
   
        function [ img ] = CTbackprojection(obj,mnidx,mxidx )
            % CTBACKPROJECTION Summary of this function goes here
            %   Detailed explanation goes here
            img = zeros(obj.nPx, obj.nPx, mxidx-mnidx+1, 'single');
            [x,y] = meshgrid(obj.x,obj.y);
            if obj.gpu == 1
                x = gpuArray(x);
                y = gpuArray(y);
            end

            for i = 1:obj.nProj
                tic
                disp(i/obj.nProj*100);
                img = img + obj.backprojection(obj.filteredProj(:,:,i),i,x,y,mnidx,mxidx);
                toc
            end
            %img(img<0)=0;
            % figure;
            % plot(t);
            % title(strcat(num2str(mean(t)),'\pm',num2str(std(t))));

        end    
        
        function vol = backprojection(obj,proj,i,xx,yy,mnidx,mxidx)
            t = obj.theta;
            if obj.dmdy>=0
               beta = (360-t(i))/360*2*pi;
            else
               beta = t(i)/360*2*pi;
            end
            vol = gpuArray(zeros(obj.nPx,obj.nPx,mxidx-mnidx+1,'single'));
            
            if obj.gpu==1
                xxS = gpuArray(xx-obj.AoRcentre);
            else
                xxS = xx-obj.AoRcentre;
            end

            u = obj.D.*(xxS*cos(beta)+yy*sin(beta)-obj.opticCentre(1)+obj.AoRposition(i))...
                   ./(xxS*sin(beta)-yy*cos(beta)+obj.D)-obj.AoRcentre+obj.opticCentre(1);

            U = (xxS*sin(beta)-yy*cos(beta)+obj.D)./obj.D;

            if obj.gpu == 1
                u = gpuArray(u);
                proj = gpuArray(proj.');
                ys = gpuArray(yy(1,1));
            if strcmp(obj.type,'f')
                if obj.dmdy>=0
                    dI = gpuArray((1+(xxS*sin(beta)-yy*cos(beta))*(-1*obj.dmdP)));
                else
                    dI = gpuArray((1+(xxS*sin(beta)-yy*cos(beta))*obj.dmdP).^2);
                end
            end
            else
                if strcmp(obj.type,'f')
                    if obj.dmdy>=0
                        dI = (1+(x*sin(beta)-y*cos(beta))*(-1*obj.dmdP)).^2;
                    else
                        dI = (1+(x*sin(beta)-y*cos(beta))*obj.dmdP).^2;
                    end
                end
            end

            for iz = mnidx:mxidx          
                if obj.gpu == 1
                    v2 = (1./U.*(-obj.nPx/2-obj.opticCentre(2)+iz))+obj.opticCentre(2);
                    
                    u2 = u.*cos(-obj.AoRangle/180*pi)-v2.*sin(-obj.AoRangle/180*pi);
                    v2 = u.*sin(-obj.AoRangle/180*pi)+v2.*cos(-obj.AoRangle/180*pi);
                    
                    u2 = u2-xxS(1,1)+1;
                    v2 = v2-ys;
                    
                    z = interp2(proj,u2,v2,obj.interptype);
                    if strcmp(obj.type,'f')
                        vol(:,:,iz-mnidx+1) = gather((dI.*z.*obj.D./U.^2).');
                    else
                        vol(:,:,iz-mnidx+1) = (z.*obj.D./U.^2).';
                        vol(:,:,iz-mnidx+1) = (z.*obj.D./U.^2).';
                    end
                else
                    v = single(frac.*obj.z(iz)+obj.opticCentre(2));
                    vol(:,:,iz-mnidx+1) = obj.D^2.*frac./(obj.D^2+v.^2).*interp2(proj',u-obj.x(1)+1,v-obj.y(1)+1,obj.interptype);
                end 


            end
            %vol(isnan(vol))=0;

            return

        end
        
        function slice = CTbackprojectionSlice(obj,sliceNo,xxS,yyS,t)
            
            slice = gpuArray(zeros(obj.nPx,obj.nPx));
            for i = 1:obj.nProj
                if obj.gpu==1
                    projI = gpuArray(obj.filteredProj(:,:,i)).';
                    projI = gpuArray(obj.projections(:,:,i)).';
                else
                    projI = (obj.filteredProj(:,:,i)).';
                    projI = (obj.projections(:,:,i)).';
                end
                
                if obj.dmdy>=0
                   beta = (360-t(i))/360*2*pi;
                else
                   beta = t(i)/360*2*pi;
                end

                op = obj.opticCentre;
                gamma = -obj.AoRangle/180*pi;
                op(1) = op(1)*cos(gamma)+op(2)*sin(gamma);
                op(2) = op(2)*cos(gamma)-op(1)*sin(gamma);
                
               
               u = obj.D.*(xxS*cos(beta)+yyS*sin(beta)-op(1)./obj.SubSampFc+obj.AoRposition(i))...
                   ./(xxS*sin(beta)-yyS*cos(beta)+obj.D)-obj.AoRcentreX/obj.SubSampFc+op(1)./obj.SubSampFc;
               
               U = (xxS*sin(beta)-yyS*cos(beta)+obj.D)./obj.D;
               
               v = (1./U.*(-obj.nPx/2-op(2)./obj.SubSampFc+sliceNo))+op(2)./obj.SubSampFc;                       
                  
                if strcmp(obj.type,'f')
                    if obj.dmdy>=0
                        dI = (1+((xxS-obj.AoRcentreX/obj.SubSampFc)*sin(beta)-yyS*cos(beta))*(-1*obj.dmdP));
                    else
                        dI = (1+((xxS-obj.AoRcentreX/obj.SubSampFc)*sin(beta)*sin(beta)-yyS*cos(beta))*obj.dmdP);
                    end
                end

                %subplot(1,2,1); imagesc(u2-u); subplot(1,2,2); imagesc(v2-v); drawnow;
                
                
                xxR = xxS.*cos(gamma)-yyS.*sin(gamma);
                yyR = yyS.*cos(gamma)+xxS.*sin(gamma);
            
                xxR = xxR-xxS(1,1)+1;
                yyR = yyR-yyS(1,1)+1;
                
                projI = interp2(projI,xxR,yyR);

    
                u2 = u-xxS(1,1)+1;
                v2 = v-yyS(1,1)+1;
                
               % u2 = u-xxS(1,1)+1;
                %v2 = v-yyS(1,1)+1;
                %subplot(2,2,1); imagesc(u2); colorbar; subplot(2,2,2); imagesc(v2); colorbar;
                
                z = interp2(projI,u2,v2,obj.interptype);
                
                if strcmp(obj.type,'f')
                    plane = dI.*z.*obj.D./U.^2;
                else
                    plane = z.*obj.D./U.^2;
                end
                %subplot(2,2,3); imagesc(projI); 
                plane(isnan(plane))=0;
                slice = slice + plane;
               % subplot(2,2,4); imagesc(plane); title(num2str(i)); drawnow;
            end
        end
        
        function slice = CTbackprojectionSlicev2(obj,sliceNo,xxS,yyS,t)   
            slice = gpuArray(zeros(obj.nPx,obj.nPx));
            op(1) = obj.opticCentre(1)*cos(obj.AoRangle/180*pi)+obj.opticCentre(2)*sin(obj.AoRangle/180*pi);
            op(2) = obj.opticCentre(2)*cos(obj.AoRangle/180*pi)-obj.opticCentre(1)*sin(obj.AoRangle/180*pi);
            xmotion = getXmotion(obj,0.033);
            for i = 1:obj.nProj
                if obj.gpu==1
                    projI = gpuArray(obj.filteredProj(:,:,i)).';
                    %projI = gpuArray(obj.projections(:,:,i)).';
                else
                    projI = (obj.filteredProj(:,:,i)).';
                    %projI = (obj.projections(:,:,i)).';
                end
                
                if obj.dmdy>=0
                   beta = (360-t(i))/360*2*pi;
                else
                   beta = t(i)/360*2*pi;
                end                
                
                s = obj.stagemotion(i)/obj.pxSz*obj.mtz0*obj.SubSampFc;
                gamma = obj.AoRangle/180*pi;
                opShift(1) = op(1)-s*(1-cos(gamma));
                opShift(2) = op(2)+s*sin(gamma);
                
                xP = 1.3325.*circshift(obj.piezomotion-mean(obj.piezomotion),[0,obj.nProj/4-1])/obj.pxSz*obj.mtz0;
                %[~,xP,~] = fitStageMotion(obj.stagemotion,obj.stepsize);
                %xP = xP/obj.pxSz*obj.mtz0;
                %shift = xP-obj.stagemotion/obj.pxSz*20.2080+6.60+obj.fullpath;
                
                shift = xP-obj.stagemotion/obj.pxSz*obj.mtz0;
                
               % shift = 0.067*sin(obj.theta/180*pi+2.595-2*pi/3+pi/2);
                %shift = (xmotion-obj.stagemotion)/obj.pxSz*obj.mtz0;
                %shift = circshift(shift,[0,200]);
                shift = shift(i);
               % shift = 0;
                %s2 = 500*sin(obj.theta/180*pi);
                %shift = obj.peaksXsub(i);
                %shift = obj.fullpath(i);
                D = obj.D;
                u = D.*(xxS*cos(beta)+yyS*sin(beta)-opShift(1)./obj.SubSampFc+obj.AoRposition(i)/obj.SubSampFc+shift)...
                    ./(xxS*sin(beta)-yyS*cos(beta)+D)+opShift(1)./obj.SubSampFc-obj.AoRcentreX/obj.SubSampFc;
                U = (xxS*sin(beta)-yyS*cos(beta)+D)./D; 
                v = (1./U.*(-obj.nPx/2-opShift(2)./obj.SubSampFc+sliceNo))+opShift(2)./obj.SubSampFc; 
             
                  
                if strcmp(obj.type,'f') 
                    if obj.dmdy>=0
                        dI = (1+((xxS-obj.AoRcentreX/obj.SubSampFc)*sin(beta)-yyS*cos(beta))*(-1*obj.dmdP));
                    else
                        dI = (1+((xxS-obj.AoRcentreX/obj.SubSampFc)*sin(beta)*sin(beta)-yyS*cos(beta))*obj.dmdP);
                    end
                end

                
                
                
                
                %% scale proj from piezo and rotate
                %scFc = 1./(1-obj.dmtdy*(1.3330.*obj.piezomotion(i)));
%               
                urot = u.*cos(gamma)-v.*sin(gamma)-s*(1-cos(gamma));
                vrot = v.*cos(gamma)+u.*sin(gamma)+s*sin(gamma);
                
                %urot = (urot-obj.opticCentre(1)/obj.SubSampFc).*scFc+obj.opticCentre(1)/obj.SubSampFc;
                %vrot = (vrot-obj.opticCentre(2)/obj.SubSampFc).*scFc+obj.opticCentre(2)/obj.SubSampFc;
                
                %%
                
                u2 = urot-xxS(1,1)+1;
                v2 = vrot-yyS(1,1);         
                
                z = interp2(projI,u2,v2,obj.interptype);
                
                %Up = ((xxS.*cos(gamma)-yyS.*sin(gamma)).*sin(beta)-(yyS.*cos(gamma)+xxS.*sin(gamma))*cos(beta)+obj.D)./obj.D;
                
                %subplot(2,2,1); imagesc(u2); colorbar; subplot(2,2,2); imagesc(v2); colorbar;
                
                if strcmp(obj.type,'f')
                    plane = dI.*z.*obj.D./U.^2;
                else
                    plane = z.*obj.D./U.^2;
                end
                %subplot(2,2,3); imagesc(projI); 
                plane(isnan(plane))=0;
                slice = slice + plane;
                %subplot(2,2,4); imagesc(plane); title(num2str(i)); drawnow; 
            end
        end
            
        function sino_out = shiftSinogram(obj,sino)
            for i=1:size(sino,2)
                s = obj.AoRposition(i);
                sino_out(:,i) = abs(rem(s,1)).*circshift(sino(:,i),[-(obj.AoRposition(i)-rem(s,1)),0])+...
                    (1-abs(rem(s,1))).*circshift(sino(:,i),[-(obj.AoRposition(i)-rem(s,1)+1),0]);
            end
        end
        
        function coordsOut = addBeadAtOptCent(obj,coords,w)
            op = obj.RopticCentre;
            newBead = [(2*rand(1)-1)*w/2/sqrt(2);(2*rand(1)-1)*w/2/sqrt(2); op(2)];
            coordsOut = horzcat(coords,newBead);
        end
        


    end

    methods (Static)
        function coords = getInitialObject(nBeads,obj)
            DOF = 1.3325*(500e-6/(0.5^2)+6.5e-3/(50*0.5));
            r1 = repmat(DOF,1,8); r2 = repmat(DOF/2,1,8); r3=0;
            phi = 0:pi/4:(2*pi-pi/4);
            r = [r1,r3,r2];
            phi = [phi,0,phi];
            x = r.*cos(phi);
            y = r.*sin(phi);
            z = linspace(min(obj.rx)*.8,max(obj.rx)*.8,nBeads);
            for i=1:nBeads
%                 if i==1
%                     coords(1,i) = obj.fullpath(1);
%                     coords(3,i) = 0;
%                     coords(2,i) = obj.piezomotion(1);
%                     coords(1,i) = coords(1,i)*cos(obj.AoRangle/180*pi)-coords(3,i)*sin(obj.AoRangle/180*pi);
%                     coords(3,i) = coords(3,i)*cos(obj.AoRangle/180*pi)+coords(1,i)*sin(obj.AoRangle/180*pi);
%                 elseif i==2
%                     coords(1,i) = obj.fullpath(1);
%                     coords(3,i) = 0+0.01;
%                     coords(2,i) = obj.piezomotion(1)+0.01;
%                     coords(1,i) = coords(1,i)*cos(obj.AoRangle/180*pi)-coords(3,i)*sin(obj.AoRangle/180*pi);
%                     coords(3,i) = coords(3,i)*cos(obj.AoRangle/180*pi)+coords(1,i)*sin(obj.AoRangle/180*pi);
%                 else
%                     coords(2,i) = (2*rand(1)-1)*max(obj.rx)-obj.piezomotion(1)*6.5e-3/obj.mtz0*.66;
%                     coords(1,i) = (2*rand(1)-1)*max(obj.rx)*.66;
%                     coords(3,i) = (2*rand(1)-1)*max(obj.rx)*.66;
%                 end          
                coords(1,i) = x(i);
                coords(3,i) = z(i);
                coords(2,i) = 1.3325*(obj.piezomotion(1)-mean(obj.piezomotion))+y(i);
                %creates beads randomly (i.e can overlap);
            end
        end

        function out = scalex(x,y,obj)
            op = obj.RopticCentre;
            dm = 1./(1-obj.dmdy.*y);
            out = (x-op(1)).*dm+op(1);
        end

        function out = scalez(z,y,obj)
            op = obj.RopticCentre;
            dm = 1./(1-obj.dmdy.*y);
            out = (z-op(2)).*dm+op(2);
        end
        
        function [xs,zs] = scaleCoords(x,y,z,obj)
            op = obj.RopticCentre;
            r =sqrt((x-op(1)).^2+(z-op(2)).^2);
            theta = atan2((z-op(2)),(x-op(1)));
            dm = 1./(1-obj.dmdy.*y);
            dr = r.*(dm-1);
            xs = x+dr.*cos(theta);
            zs = z+dr.*sin(theta);
        end

        function [x,y,z] = coordsMagChange(x,y,z,obj)
                %find depth position and calculates relative change in magnification
                %and then scales the x and z positions respectively.
                dm = 1./(1-obj.dmdy.*y);
                op = obj.RopticCentre;
                x = (x-op(1)).*dm+op(1);
                z = (z-op(2)).*dm+op(2);
        end
        
        function out = scalexPixels(x,y,obj)
            op = obj.opticCentre;
            dm = 1./(1-obj.dmdP.*y);
            out = (x-op(1)).*dm+op(1);
        end
        
        function out = scalezPixels(z,y,obj)
            op = obj.opticCentre;
            dm = 1./(1-obj.dmdP.*y);
            out = (z-op(1)).*dm+op(1);
        end
        
        function [x,y,z] = coordsMagChangeOpp(x,y,z,obj)
                %find depth position and calculates relative change in magnification
                %and then scales the x and z positions respectively.
                dm = (1+obj.dmdy.*y);
                op = obj.RopticCentre;
                x = (x-op(1)).*dm+op(1);
                z = (z-op(2)).*dm+op(2);
        end

        function coordsOut = rotateSimp(coords,theta)
        %rotates coords about AoR centre. 0 is middle of image. Need to
        %convert beta to radians
        theta = theta/180*pi;
        coordsOut(1,:) = (coords(1,:)).*cos(theta)-coords(2,:).*sin(theta);
        coordsOut(2,:) = (coords(1,:)).*sin(theta)+coords(2,:).*cos(theta);
        coordsOut(3,:) = coords(3,:);
        end
        
        %moves coord array, translation matrix moves the point to (0,0,0)
        function coordsOut = T(coords,point)
            coordsOut(1,:) = coords(1,:) - point(1);
            coordsOut(2,:) = coords(2,:) - point(2);
            coordsOut(3,:) = coords(3,:) - point(3);      
        end
        
        %translate xx,yy meshgrid by dx, dy
        function [xxOut,yyOut] = Tm(xx,yy,dx,dy)
           xxOut = xx + dx;
           yyOut = yy + dy;
        end
        
        %Rotation of coords around unit vector (u,v,w) (u^2+v^2+w^2)=1 by
        %angle theta
        function coordsOut = rotateComp(coords,u,v,w,theta)
            theta = theta/180*pi;
            x = coords(1,:);
            y = coords(2,:);
            z = coords(3,:);
            coordsOut(1,:) = u*(u.*x+v.*y+w.*z)*(1-cos(theta))+x*cos(theta)+(-w.*y+v.*z)*sin(theta);
            coordsOut(2,:) = v*(u.*x+v.*y+w.*z)*(1-cos(theta))+y*cos(theta)+(w.*x-u.*z)*sin(theta);
            coordsOut(3,:) = w*(u.*x+v.*y+w.*z)*(1-cos(theta))+z*cos(theta)+(-v.*x+u.*y)*sin(theta);
        end
        
        function [h, nn] = ramp_flat(n)
            nn = [-(n/2):(n/2-1)]';
            h = zeros(size(nn),'single');
            h(n/2+1) = 1 / 8;
            odd = mod(nn,2) == 1;
            h(odd) = -0.5 ./ (pi * nn(odd)).^2;
        end

        function [filt] = Filter(filter, kernel, order, d)

        f_kernel = abs(fft(kernel))*2;
        filt = f_kernel(1:order/2+1)';
        w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

        switch lower(filter)
            case 'ram-lak'
                % Do nothing
            case 'shepp-logan'
                % be careful not to divide by 0:
                filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
            case 'cosine'
                filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
            case 'hamming'  
                filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
            case 'hann'
                filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
            otherwise
                filter
                error('Invalid filter selected.');
        end

        filt(w>pi*d) = 0;                      % Crop the frequency response
        filt = [filt , filt(end-1:-1:2)];    % Symmetry of the filter
        return
        end
        
        function [rotAngle,shiftAm] = findRotShift(mip,varargin)

if nargin<2
k=1; angles = -5:0.5:5; shifts = -100:25:100;

for i=shifts
    l=1;
    for j=angles
        s(k,l,:) = rotation(mip,i,j);
        l=l+1;
    end
    k=k+1;
    disp(i);
end; 


figure;
for i=1:4
    subplot(2,2,i);
    imagesc(angles,shifts,s(:,:,i));
end

figure; imagesc(angles,shifts,std(s,0,3));

dif = std(s,[],3);
[r,c] = find(dif==min(dif(:)));
rotAngle = angles(c); shiftAm = shifts(r);
sprintf('Initial guess at shift=%f,rot=%f,sse=%f',rotAngle,shiftAm,min(dif(:)))
else
     shiftAm = varargin{1}; rotAngle = varargin{2};
end

start_guesses = [shiftAm,rotAngle];
model = @fun;
est = fminsearch(model,start_guesses);

rotAngle = est(2); shiftAm = est(1);

    function sse = fun(params)
        A = params(1);
        B = params(2);
        if or(or(A<-100,A>100),or(B<-5,B>5))
            sse=1;
        else
        out = rotation(mip,A,B);
        sse = std(out,0,2);
        sprintf('shift=%f,rot=%f,sse=%f',A,B,sse)
        end
    end
        end  
    end
end

