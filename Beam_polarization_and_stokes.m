clc, clear, close all

% Physical parameters
w0 = 1e-3;                                     % Beam width
m = 0;                                         % Topological charge
N = 2^8;                                       % Numerical window
E = @(r,th) (r/w0).^(abs(m)).*exp(-r.^2/w0^2).*exp(1i*m*th); % Electric field (Laguerre-Gauss)

% Mesh's space generation
[x,y] = meshgrid(linspace(-2.5*w0,2.5*w0,N));  % Cartesian coordinates
[th,r] = cart2pol(x,y);                        % Polar coordinates

% Physical magnitudes
intensity = abs(E(r,th).^2);                   % Electric field's intensity
phase = angle(E(r,th));                        % Electric field's phase

% Intensity and phase graphing
figure, surf(x,y,intensity,"EdgeColor","none"), axis([-1 1 -1 1]*2.5*w0)
xlabel("x"), ylabel("y"), zlabel("E"), view(2)
title("Beam's intensity, m = "+num2str(m)), colormap gray
figure, surf(x,y,phase,"EdgeColor","none"), axis([-1 1 -1 1]*2.5*w0)
xlabel("x"), ylabel("y"), zlabel("arg(E)"), view(2)
title("Beam's phase, m = "+num2str(m)), colormap gray

%% BEAM POLARIZATION
% vertical (V), horizontal (H), diagonal(D), antidiagonal (A),
% right-circular (R), left-circular (L)
polarization = "V";
[Ex,Ey] = polarize(E(r,th),polarization);

% Beam's vector field graphing
if (polarization == "L") || (polarization == "R")
    w = 0.1;
    figure
    for t = 0:2*pi*20
        quiver(x,y,real(Ex*exp(-1i*w*t)),real(Ey*exp(-1i*w*t)),"w")
        xlabel("x"), ylabel("y"), title("Polarized beam")
        axis square, set(gca,"Color","k"), axis([-1 1 -1 1]*2*w0) 
        drawnow
    end
else
    figure
    quiver(x,y,real(Ex),real(Ey),"w")
    xlabel("x"), ylabel("y"), title("Polarized beam")
    axis square, set(gca,"Color","k"), axis([-1 1 -1 1]*2*w0) 
end



% Normalized stokes parameters
I0 = I(Ex,Ey,0,0) + I(Ex,Ey,pi/2,0)
S0 = (I(Ex,Ey,0,0) + I(Ex,Ey,pi/2,0))/I0
S1 = (I(Ex,Ey,0,0) - I(Ex,Ey,pi/2,0))/I0
S2 = (I(Ex,Ey,0,pi/4) - I(Ex,Ey,0,-pi/4))/I0
S3 = (I(Ex,Ey,pi/2,pi/4) + I(Ex,Ey,pi/2,-pi/4))/I0

% Polarization grade
P = sqrt(S1^2 + S2^2 + S3^2)/S0


function [Ex,Ey] = polarize(E,Polarization)
    switch Polarization
        case "H"
            p = [1;0];              % Horizontal polarization
        case "V"
            p = [0;1];              % Vertical polarization
        case "D"
            p = 1/sqrt(2)*[1;1];    % Diagonal polarization
        case "A"
            p = 1/sqrt(2)*[1;-1];   % Anti-diagonal polarization
        case "R"
            p = 1/sqrt(2)*[1;-1i];  % Right-circular polarization
        case "L"
            p = 1/sqrt(2)*[1;1i];   % Left-circular polarization
    end
    [Ex,Ey] = deal(E*p(1),E*p(2));  % Electric field's components polarization
end

function [Ex,Ey] = retarder(Ex,Ey,phi)  % Function representing a retarder element
    Ex = Ex;
    Ey = Ey*exp(1i*phi);
end

function [Ex,Ey] = polarizator(Ex,Ey,angle) % Function representing a polarizator element
    Ex = Ex*sin(angle).^2 + Ey*sin(angle)*cos(angle);
    Ey = Ex*sin(angle)*cos(angle) + Ey*cos(angle).^2;    
end

function Intensity = I(Ex,Ey,theta, phi)
    [Ex,Ey] = retarder(Ex,Ey,phi);              % Beam retardation
    [Ex,Ey] = polarizator(Ex,Ey,theta);         % Beam polarization
    magEx = sum(Ex,"all")*conj(sum(Ex,"all"));
    magEy = sum(Ey,"all")*conj(sum(Ey,"all"));
    Intensity = magEx + magEy;                  % Beam's intensity
end