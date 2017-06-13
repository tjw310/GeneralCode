% zero = coneRecon([0,0],0);
% om= coneRecon([-50,0],0);
% op= coneRecon([50,0],0);
% rm= coneRecon([0,0],-20);
% rp= coneRecon([0,0],20);
% oprp= coneRecon([50,0],20);
% oprm= coneRecon([50,0],-20);
% omrp= coneRecon([-50,0],20);
% omrm= coneRecon([-50,0],-20);
% 
% sz = 3e-3;

%%
zero.dmdy=-10;
om.dmdy=-10;
op.dmdy=-10;
rm.dmdy=-10;
rp.dmdy=-10;
oprp.dmdy=-10;
oprm.dmdy=-10;
omrp.dmdy=-10;
omrm.dmdy=-10;
m.dmdy=-10;
m.AoRmotion = -10*sin(m.theta/180*pi);

zero.testObject(2,sz,zero.coordsInit);
om.testObject(2,sz,zero.coordsInit);
%%
op.testObject(2,sz,zero.coordsInit);
rm.testObject(2,sz,zero.coordsInit);
rp.testObject(2,sz,zero.coordsInit);
m.testObject(2,sz,zero.coordsInit);
oprp.testObject(2,sz,zero.coordsInit);
oprm.testObject(2,sz,zero.coordsInit);
omrp.testObject(2,sz,zero.coordsInit);
omrm.testObject(2,sz,zero.coordsInit);
%%
%zero.fanBeamrecon(squeeze(zero.projections(:,128,:)));
om.fanBeamrecon(squeeze(om.projections(:,128,:)));
op.fanBeamrecon(squeeze(op.projections(:,128,:)));
rm.fanBeamrecon(squeeze(rm.projections(:,128,:)));
rp.fanBeamrecon(squeeze(rp.projections(:,128,:)));
oprp.fanBeamrecon(squeeze(oprp.projections(:,128,:)));
oprm.fanBeamrecon(squeeze(oprm.projections(:,128,:)));
omrp.fanBeamrecon(squeeze(omrp.projections(:,128,:)));
omrm.fanBeamrecon(squeeze(omrm.projections(:,128,:)));