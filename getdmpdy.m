function [dmpdy,varargout] = getdmpdy(deltaM,yp2,yp1,varargin)
%gets the gradient of piezo mag change. deltaM = x2/x1 where x2 is the x
%location for yp2, and x1 is the x location for yp1. yp2 and yp1 are the
%piezo y-locations (and the actual y-locations of the beads). This is based
%on keeping the ETL constant and changing the piezo position.

%varagin arguments: errDeltaM, errYp2, errYp1
%varargout arguments: error in dmpdy

dmpdy = (deltaM-1)./(deltaM.*yp2-yp1);

if and(nargin>3,nargin==6)
    errDeltaM = varargin{1};
    errYp2 = varargin{2};
    errYp1 = varargin{3};
    
    errordmpdy = 1./(deltaM*yp2-yp1).^2*((yp2-yp1).^2.*errDeltaM.^2+(deltaM-1).^2.*(errYp1.^2+deltaM.^2.*errYp2.^2))^.5;
    varargout{1} = errordmpdy;
elseif and(nargin>3,nargin~=6)
    error('enter correct number of arguments');
end
    
end