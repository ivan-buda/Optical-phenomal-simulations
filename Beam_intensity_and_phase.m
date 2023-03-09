clc, clear, close all

% Physical parameters
w0 = 1e-3;                  % Beam width
m = 0;                      % Topological charge
E = @(r,th) (r/w0).^(abs(m)).*exp(-r.^2/w0^2).*exp(1i*m*th); % Electric field

% Mesh's space generation
dr = w0/100; dth = 0.01*2*pi;
rr = 0:dr:2.5*w0; thth = 0:dth:2*pi+dth;
[r,th] = meshgrid(rr,thth); % Cartesian coordinates
[x,y] = pol2cart(th,r);     % Polar coordinates

% Physical magnitudes
I = abs(E(r,th).^2);        % Electric field's intensity
phase = angle(E(r,th));     % Electric field's phase

% Intensity and phase graphing
figure, surf(x,y,I,"EdgeColor","none"), axis([-1 1 -1 1]*2.5*w0)
xlabel("x"), ylabel("y"), zlabel("E"), view(2)
title("Electric field's intensity, m = "+num2str(m)), colormap gray
figure, surf(x,y,phase,"EdgeColor","none"), axis([-1 1 -1 1]*2.5*w0)
xlabel("x"), ylabel("y"), zlabel("arg(E)"), view(2)
title("Electric field's phase, m = "+num2str(m)), colormap gray

% Initial intensity
I0 = sum(sum(abs(E(r,th)).^2))*dr*dth;

% MAlus law
theta = 0:0.1:2*pi;
I = @(theta) I0*(cos(theta)).^2;
figure, plot(theta,I(theta)), xlim([0 max(theta)])
title("Malus law"), xlabel("\theta [\pi rad]"), ylabel("I(\theta)")