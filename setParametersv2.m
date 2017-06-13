function [param] = setParametersv2()
%General Parameters
%scan range
param.SR = 1.1; %mm (ful range)
%n.projections
param.nProj = 100;
%angles covered
param.nAngles = [0,360];
%angles
param.theta = param.nAngles(1):(param.nAngles(2)-param.nAngles(1))/param.nProj:param.nAngles(2)-(param.nAngles(2)-param.nAngles(1))/param.nProj;
%object space size in pixels;
param.nPx = 256;
%Field of view
param.FOV = 128; %in pixels
%subsampling factor;
param.SubSampFc = param.FOV/param.nPx;
%effeective pixel size
param.pxSz = 6.5e-3*param.SubSampFc; %mm

%% parameteres from experimental data used to calcuate the cone height. 
%object rotates around the z-axis. y-axis is depth
%gradient of static magnification of fixed plane.
param.dmtdy = 0.2220;
%transverse magnification at z=0
param.mtz0 = 20.3063;
%real size
param.w = param.nPx*param.pxSz/param.mtz0; %real size in mm of object space
param.SR = param.w; %set object space to square
%relative magnification gradient of object
param.dmdy = -0.145;
param.dmdy = -3;
param.dmdP = param.dmdy.*param.pxSz./param.mtz0; %relative mag change per z pixel (z goes about 0)
%cone height effective source-detector distance
param.Dm = -1/param.dmdy;
param.D = -1/param.dmdP;
% half cone angle
param.halfConeAngle = atan2(param.w/2,param.Dm)/pi*180; %in degrees

%%
%optical centre in real size(centre of rotation in middle of volume)
param.opticCentreP = [-75,-50]; %in pixels
param.opticCentreP = [0,0]; %in pixels
param.opticCentreP = [75,50]; %in pixels
param.opticCentre = param.opticCentreP*param.pxSz./param.mtz0; %in mm
%detector cone AoR offset
param.off_u = param.opticCentreP(1); param.off_v = param.opticCentreP(2); % detector rotation shift (pixels)

%distance from source to AoR;
param.R = sqrt(param.opticCentreP(1)^2+param.D^2); %pixels

%Object coordinates
param.x = -param.nPx/2:param.nPx/2-1; % number of voxels
param.y = -param.nPx/2:param.nPx/2-1;
param.z = -param.nPx/2:param.nPx/2-1;

%real scales
param.rx = param.x.*param.pxSz./param.mtz0;
param.ry = param.y.*param.pxSz./param.mtz0;
param.rz = param.z.*param.pxSz./param.mtz0;

%The real detector panel pixel density (number of pixels)
param.u =  param.x - param.off_u;		% number of pixels
param.v =  param.z;

param.parker = 0; % data with 360 deg -> param.parker = 0 , data less than 360 deg -> param.parker=1 

param.interptype = 'linear'; % 'linear', 'nearest'

% % % % % % Confirm your parameters % % % % % % %
% Only for Matlab version above 2013b with parallel computing toolbox: Speed dependent on your GPU
% You don't need to install CUDA, but install latest graphics driver.
% only Nvidia GPU cards can use this. otherwise please "param.gpu=0"
% This option is semi-GPU code using built-in Matlab GPU functions: several times faster than CPU
param.gpu = 1;

param.filter = 'ram-lak';
end


