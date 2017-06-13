function g = czt_TW2(xz,dx,dz,m,varargin)
%CZT  Chirp z-transform, given input field x with associated scale increment dx, f, output size of m
% Version 2.0 allows the setting of 4 frequency limits, 2 for each
% dimension. Only 2D functions allowed.

uMax = 1/(2*dx);
vMax = 1/(2*dz);

if nargin>5
    if nargin==6
        fMin = varargin{1};
        fMax = varargin{2};
        gMin = varargin{1};
        gMax = varargin{2};
    elseif nargin==8
        fMin = varargin{1};
        fMax = varargin{2};
        gMin = varargin{3};
        gMax = varargin{4};
    else
        error('Must have 5,6 or 8 arguments');
    end
        
elseif nargin==5
    fMin = -varargin{1};
    fMax = varargin{1};
    gMin = -varargin{1};
    gMax = varargin{1};
else
    error('Must have 5 arguments');
end
    
w2 = exp(-1i*2*pi*(fMax-fMin)/(2*m*uMax));
a2 = exp(1i*2*pi*fMin/(2*uMax));       

w = exp(-1i*2*pi*(gMax-gMin)/(2*m*vMax));
a = exp(1i*2*pi*gMin/(2*vMax));   

g = czt_JC(czt_JC(xz,m,w2,a2).',m,w,a).';  


end
    
    



