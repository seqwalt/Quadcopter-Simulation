% --- Nonlinear dynamics of the quadcopter --- %

%function dydt = quad_dyn(t,y,pd,xref,umax,refmax)
function dydt = quad_nonlinear(t,y,K,ref)

g = -9.81; % m/s^2
l = .15;   % meters
m = .4;    % kg
Ix = 0.25*m*l^2;
Iy = Ix;
Iz = .5*m*l^2;

% State variables
phi = y(1);
th = y(2);
psi = y(3);
p = y(4);
q = y(5);
r = y(6);
u = y(7);
v = y(8);
w = y(9);
X = y(10)-ref(1);
Y = y(11);
Z = y(12)-ref(2);

% Wind disturbances
fwx = 0;
fwy = 0;
fwz = 0;
tauwx = 0;
tauwy = 0;
tauwz = 0;

% Controller
STATE = [phi; th; psi; p; q; r; u; v; w; X; Y; Z];
U = -K*STATE;
ft = m*g - U(1);
taux = U(2);
tauy = U(3);
tauz = U(4);

% Equations of motion
dphi = p + r*cos(phi)*tan(th) + q*sin(phi)*tan(th);
dth = q*cos(th) - r*sin(phi);
dpsi = [r*cos(phi)/cos(th)] + [q*sin(phi)/cos(th)];
dp = [(Iy - Iz)*r*q/Ix] + [(taux + tauwx)/Ix];
dq = [(Iz - Ix)*p*r/Iy] + [(tauy + tauwy)/Iy];
dr = [(Ix - Iy)*p*q/Iz] + [(tauz + tauwz)/Iz];
du = r*v - q*w - g*sin(th) + fwx/m;
dv = p*w - r*u + g*sin(phi)*cos(th) + fwy/m;
dw = q*u - p*v + g*cos(th)*cos(phi) + (fwz-ft)/m;
dx = w*[sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(th)] - ...
    v*[cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(th)] + u*cos(psi)*cos(th);
dy = v*[cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(th)] - ...
    w*[cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(th)] + u*cos(th)*sin(psi);
dz = w*cos(phi)*cos(th) - u*sin(th) + v*cos(th)*sin(phi);

dydt = [dphi; dth; dpsi; dp; dq; dr; du; dv; dw; dx; dy; dz];

% Ground detection
if y(12) < 0
    dydt = zeros(12,1);
end