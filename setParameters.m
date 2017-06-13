function [param] = setParameters()
%General Parameters
param.pxSz = 6.5e-3*6.4; %mm
%scan range
param.SR = 1.1; %mm (ful range)
%n.projections
param.nProj = 100;
%angles covered
param.nAngles = [0,360];
%angles
param.theta = param.nAngles(1):(param.nAngles(2)-param.nAngles(1))/param.nProj:param.nAngles(2)-(param.nAngles(2)-param.nAngles(1))/param.nProj;
%object size;
param.nPx = 128;


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
%cone height effective source-detector distance
param.DSD = param.SR/2-1/param.dmdy;

%%
%optical centre in real size(centre of rotation in middle of volume)
%param.opticCentre = [50,-50].*param.pxSz/param.mtz0;
param.opticCentre = [75,-50];
%detector cone AoR offset
param.off_u = param.opticCentre(1); param.off_v = param.opticCentre(2); % detector rotation shift (real size)
%param.off_u=0; param.off_v=0;

%Object coordinates
param.nx = param.nPx; % number of voxels
param.ny = param.nPx;
param.nz = param.nPx;

% X-ray source and detector setting
param.DSO = -1/param.dmdy;	%  X-ray source to object axis distance
%param.DSD = param.DSO;
%param.D = param.DSO;

param.sx = param.nPx*param.pxSz/param.mtz0; % mm (real size)
param.sy = param.nPx*param.pxSz/param.mtz0;
param.sz = param.nPx*param.pxSz/param.mtz0;

%The real detector panel pixel density (number of pixels)
param.nu = param.nPx;		% number of pixels
param.nv = param.nPx;

% Detector setting (real size)
param.su = param.nPx*param.pxSz/param.mtz0.*param.DSO/param.DSD;	% mm (real size)
param.sv = param.nPx*param.pxSz/param.mtz0.*param.DSO/param.DSD;  % mm

param.parker = 0; % data with 360 deg -> param.parker = 0 , data less than 360 deg -> param.parker=1 

% % % Geometry calculation % % %
param.xs = linspace(-param.sx/2,param.sx/2,param.nPx);
param.ys = linspace(-param.sy/2,param.sy/2,param.nPx);
param.zs = linspace(-param.sz/2,param.sz/2,param.nPx);

param.us = linspace(-param.su/2,param.su/2,param.nPx) + param.off_u;
param.vs = linspace(-param.sv/2,param.sv/2,param.nPx) + param.off_v;

param.du = param.pxSz/param.mtz0;
param.dv = param.pxSz/param.mtz0;

param.interptype = 'linear'; % 'linear', 'nearest'

% % % % % % Confirm your parameters % % % % % % %
% Only for Matlab version above 2013b with parallel computing toolbox: Speed dependent on your GPU
% You don't need to install CUDA, but install latest graphics driver.
% only Nvidia GPU cards can use this. otherwise please "param.gpu=0"
% This option is semi-GPU code using built-in Matlab GPU functions: several times faster than CPU
param.gpu = 1;

param.filter = 'ram-lak';
end


