% --- Nonlinear dynamics of the quadcopter --- %

%function dydt = quad_dyn(t,y,pd,xref,umax,refmax)
function dydt = quad_linear(t,y,A,B,K,ref)

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

% Controller
STATE = [phi; th; psi; p; q; r; u; v; w; X; Y; Z];
U = -K*STATE;

% Equations of motion
dydt = A*STATE+B*U;

% Ground detection
if y(12) < 0
    dydt = zeros(12,1);
end