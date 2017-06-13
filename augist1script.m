% 1st august 2016 script

ft = 180;
fo = 9;
beta = 3;
fd = -75;

Alpha = 7:0.01:11;
Gamma = 2:0.01:10;

fe = 45:0.1:120;
I = 1./fe;

I2 = linspace(0,290,length(fe));

I3 = (I-min(I))./(max(I)-min(I))*length(fe)+1;

fe = interp1(fe,I3);
clear M
k=1;
for gamma = Gamma

    A = -ft/fd+gamma*ft/fd./fe+ft*beta/fo/fd-gamma*beta*ft./fe/fo/fd+gamma*ft/fo/fd-ft./fe+beta*ft./fe/fo-ft/fo;

    M(k,:) = A;
    k=k+1;

end

figure; surf(I2,Gamma,M,'edgecolor','none'); xlabel('ETL Current'); ylabel('Piezo Shift'); zlabel('Mag'); colorbar;