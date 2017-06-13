function runGroupRecon(path,angle,offset,motion,images)

%path = uigetdir;

p1 = fullfile(path,'ETLonly');
p2 = fullfile(path,'Piezo');

p1a = fullfile(p1,'data30mA');
p1b = fullfile(p1,'data60mA');

p2a = fullfile(p2,'data1mA');
p2b = fullfile(p2,'data30mA');
p2c = fullfile(p2,'data60mA');

images.path = p1a; runRecon(images,angle,offset,motion);
images.path = p1b; runRecon(images,angle,offset,motion);
images.path = p2a; runRecon(images,angle,offset,motion);
images.path = p2b; runRecon(images,angle,offset,motion);
images.path = p2c; runRecon(images,angle,offset,motion);

function runRecon(images,angle,offset,motion)
images.piezomotion = zeros(1,400);
images.projections = [];
images.filteredProj = [];
images.loadProjections(images.path);
images.AoRangle = angle;
images.projections = scaleProjections(images);
images.AoRangle = 0;
images.AoRmotion = motion;
images.AoRcentreX = offset;
images.reconPS(900,1615);
end

end

