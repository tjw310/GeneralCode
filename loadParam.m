%% Parameter setting  TW Edit%%

% % % % % % Confirm your parameters % % % % % % %

function param = loadParam(images,magCalib,mappedTracks)
% magCalib is the magnificationCalibration class object, defining the
% system and mappedTracks is the associated mapped bead tracks that will
% get the scan range and magnification change of the system

%SR - scan range and dmrdz is the relative magnfication change of
%beads(including the transverse mag change).
param.SR = mean(magCalib.scatterzSR(mappedTracks));
param.dmrdz = mean(magCalib.scatterZmagGrad(mappedTracks));

%dmtdz is the gradient of the fixed mag change, and mt is the magnification
%at z=0. This is from measurements of a grating.
param.mt = 20.3063;
param.dmtdz = 0.2220;

%beta is the source-detector distance. see PDF for maths.
param.beta = (param.dmrdz*param.SR/2-1)/(param.dmrdz-1/param.mt*param.dmtdz);

%F is the size of the field of view at maximum magnfication. G is the size
%of the field of view at minimum magnification.
param.F = size(images,1);
param.G = param.F*(1-param.SR/param.beta);

%pxSz - real pixel size (in mm)
param.pxSz = 6.5e-3;

param.nx = param.G; % number of voxels
param.ny = param.G;
param.nz = param.G;

param.sx = param.G*param.pxSz/(param.mt+param.dmtdz*param.SR/2); % mm (real size)
param.sy = param.G*param.pxSz/(param.mt+param.dmtdz*param.SR/2); % mm
param.sz = param.G*param.pxSz/(param.mt+param.dmtdz*param.SR/2); % mm

%The real detector panel pixel density (number of pixels)
param.nu = size(images,1);		% number of pixels
param.nv = size(images,2);

% Detector setting (real size)
param.su = 1024;	% mm (real size)
param.sv = 800;     % mm

% X-ray source and detector setting
param.DSD = param.beta;    %  Distance source to detector 
param.DSO = param.beta-param.SR/2;	%  X-ray source to object axis distance

% angle setting
param.nProj = size(images,3);
param.dir = 1;   % gantry rotating direction (clock wise/ counter clockwise)
param.dang = 360/param.nProj; % angular step size (deg)
param.deg = 0:param.dang:360-param.dang; % you can change
param.deg = param.deg*param.dir;


param.parker = 0; % data with 360 deg -> param.parker = 0 , data less than 360 deg -> param.parker=1 

% % % % % % Confirm your parameters % % % % % % %
 
% filter='ram-lak','cosine', 'hamming', 'hann' 
param.filter='ram-lak'; % high pass filter

param.dx = param.sx/param.nx; % single voxel size
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

param.off_u = 0; param.off_v = 0; % detector rotation shift (real size)

% % % Geometry calculation % % %
param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx;
param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz;

param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v;

param.interptype = 'linear'; % 'linear', 'nearest'

% % % % % % Confirm your parameters % % % % % % %
% Only for Matlab version above 2013b with parallel computing toolbox: Speed dependent on your GPU
% You don't need to install CUDA, but install latest graphics driver.
% only Nvidia GPU cards can use this. otherwise please "param.gpu=0"
% This option is semi-GPU code using built-in Matlab GPU functions: several times faster than CPU
param.gpu = 1;

end


















