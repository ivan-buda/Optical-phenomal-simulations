clc, clear, close all

%% Shoot photons with a gaussian density
Nphotons = 1000e3;      % Number of photons shooted
[xdata,ydata] = deal(randn(1,Nphotons),randn(1,Nphotons)); % Photons location

%% Create histogram of the photons' location
bins = -5:0.03:5;
[X,Y] = meshgrid(bins,bins);
XYdata = [xdata;ydata]';

I = hist3(XYdata,{bins bins});
xcut = I(end/2,:);

%% Plot results
figure
plot(xdata,ydata,"w.","MarkerSize",0.01)
title("Photons collided on screen"); 
xlabel("x")
ylabel("y")
set(gca,"Color","k")
axis([min(xdata) max(xdata) min(ydata) max(ydata)]*1.5)

figure
surf(X,Y,I)
title("Gaussian distribution of photons")
xlabel("x")
ylabel("y")
zlabel("Intensity")
view(2)
colormap hot
shading interp

figure
plot(bins,xcut,"r")
title("Gaussian distribution in the x axis")
xlabel("x")
ylabel("y")
grid on