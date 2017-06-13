classdef gaussFit
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a
        b
        c
        d
    end
    
    methods
        function obj = gaussFit(a,b,c,d)
            obj.a=a;
            obj.b=b;
            obj.c=c;
            obj.d=d;
        end
        
        function out = plot(obj,x,pxSz,varargin)
            out = obj.a + (obj.b-obj.a).*exp(-(x-obj.c).^2/(2*obj.d^2));
            if nargin>3
                plot(x.*pxSz,out);
                drawnow;
            end
        end
        
        function out = normplot(obj,x,pxSz,varargin)
            out = obj.a + (obj.b-obj.a).*exp(-(x-obj.c).^2/(2*obj.d^2));
            if nargin>3
                plot(x.*pxSz,out./max(out));
                drawnow;
            end
        end
        
        function out = grad(obj,x,pxSz,varargin)
            y = obj.plot(x);
            out = gradient(y)./pxSz;
            if nargin>3
                plot(x.*pxSz,out);
                drawnow;
            end
        end
        
        %get normalised gradient. Normalise y values before calculating
        function out = normgrad(obj,x,pxSz,varargin)
            y = obj.plot(x,pxSz);
            out = gradient(y./max(y),pxSz);
             if nargin>3
                plot(x.*pxSz,out);
                drawnow;
             end
        end
        
        function out = errornormgrad(obj,x,pxSz,varargin)
            dy = obj.normgrad(x,pxSz);
            out = 
        end
            
    end
    
end

