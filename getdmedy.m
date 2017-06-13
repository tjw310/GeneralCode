function [dmedy,varargout] = getdmedy(deltaM,yp,y2,y1,varargin)
%gets the gradient of ETL mag change. deltaM = x2/x1 where x2 is the x
%location for ye2, and x1 is the x location for ye1. ye2 and ye1 are the
%actual y-locations, and yp is the piezo location. This is based on keeping
%the piezo location fixed and changing the ETL focal positions

%varagin arguments: errDeltaM, errYp, errY2, errY1
%varargout arguments: error in dmpdy

dmedy = (deltaM-1)./(deltaM.*(y2-yp)-(y1-yp));

if and(nargin>3,nargin==8)
    errDeltaM = varargin{1};
    errYp = varargin{2};
    errY2 = varargin{3};
    errY1 = varargin{4};
    
    v = deltaM.*(y2-yp)-(y1-yp);
    
    errordmedy = 1./v.^2.*((y2-y1).^2.*errDeltaM.^2+(deltaM-1).^2.*((deltaM-1).^2.*errYp.^2+deltaM.^2.*errY2.^2+errY1.^2)).^.5;   
    varargout{1} = errordmedy;
elseif and(nargin>3,nargin~=8)
    error('enter correct number of arguments');
end

end