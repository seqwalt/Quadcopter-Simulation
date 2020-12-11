%% --- ECE 717 Final Project --- %%

g = sym('g');
m = sym('m');
Ix = sym('Ix');
Iy = sym('Iy');
Iz = sym('Iz');
syms t t1 refx refz;

A = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0,-g, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     g, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];

B = [0,    0,    0,    0;
     0,    0,    0,    0;
     0,    0,    0,    0;
     0, 1/Ix,    0,    0;
     0,    0, 1/Iy,    0;
     0,    0,    0, 1/Iz;
     0,    0,    0,    0;
     0,    0,    0,    0;
   1/m,    0,    0,    0;
     0,    0,    0,    0;
     0,    0,    0,    0;
     0,    0,    0,    0];

%% State-Space Analysis

% Stability -> Unstable
% All zeros eig vals, and D has more cols than V, implying the
% eigenvectors are lin dependent -> not stable. Because all eigenvalues of
% A are zeros, then the number of eigenvalues must equal the size of the
% null space of A for stability, but this is not true.
[V,D] = eig(A);

% Controllability -> Controllable
% rank(Ctr) = n = 12 is satisfied
Ctr = [B A*B A*A*B pow(A,3)*B pow(A,4)*B pow(A,5)*B pow(A,6)*B ...
       pow(A,7)*B pow(A,8)*B pow(A,9)*B pow(A,10)*B pow(A,11)*B];
rank_ctr = rank(Ctr);

%% Controller Design

g = -9.81; % m/s^2
l = .15;   % meters
m = .5;    % kg
Ix = 0.25*m*l^2;
Iy = Ix;
Iz = .5*m*l^2;
state_ref = [0; 0; 0; 0; 0; 0; 0; 0; 0; refx; 0; refz];
A_num = double(subs(A));
B_num = double(subs(B));

% Feedforward Control
expr = expm(A_num*t)*(B_num*B_num')*expm(A_num'*t);
Wr = int(expr,t,0,t1);
span = rref(Wr); % reduced row echelon form of Wr will show what vectors make up its span
det(Wr); % not equal to zero when t0 does not equal t1 -> Wr is invertable
eta = inv(Wr)*state_ref;
u_closed = B_num'*expm(A_num'*(t1))*eta;

% Optimal Feedback Control - LQR
Q = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % phi
     0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % th
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;  % psi
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;  % p
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;  % q
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;  % r
     0, 0, 0, 0, 0, 0,10, 0, 0, 0, 0, 0;  % u
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;  % v
     0, 0, 0, 0, 0, 0, 0, 0,10, 0, 0, 0;  % w
     0, 0, 0, 0, 0, 0, 0, 0, 0,10, 0, 0;  % x
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;  % y
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,10]; % z
 
R = [2, 0, 0, 0;  % ft
     0, 1, 0, 0;  % taux
     0, 0, 2, 0;  % tauy
     0, 0, 0, 1]; % tauz
 
[K,S,e] = lqr(A_num,B_num,Q,R);
K(abs(K) < 1e-8) = 0;

% Stability of controlled system -> Stable
% rank(V_cl) = num of eigenvalues
A_cl = A_num + B_num*K;
[V_cl,D_cl] = eig(A_cl);

%% Simulation

LABEL = {"","z (m)",["Linear Simulation","a1)"],"","",["Nonlinear Simulation","b1)"];
         "","z (m)","a2)","","","b2)";
         "","z (m)","a3)","","","b3)";
         "x (m)","z (m)","a4)","x (m)","","b4)"};

t_l = 6;   % sim time length in seconds
freq = 50; % measurement frequency in Hz
tspan = linspace(0,t_l,t_l*freq);

IC  = [0,0,-4,2;
       0,0,0,2;
       0,0,0,2;
       -2,4,0,2];
REF = [0,0;
       3,2;
       6,3.5;
       2,2];
rows = size(IC,1);
% Feed-BACK control
figure
for i = 1:rows
    % Format of the initial conditions:
    %           [phi0; th0; psi0; p0; q0; r0;       u0; v0;       w0;       x0; y0; z0]
    init_cond = [   0;   0;    0;  0;  0;  0;  IC(i,1);  0;  IC(i,2);  IC(i,3);  0; IC(i,4)];
    ref = REF(i,:); % [x_ref,z_ref]
    [t,states] = ode45(@(t,y) quad_nonlinear(t,y,K,ref),tspan,init_cond);
    [t_lin,states_lin] = ode45(@(t,y) quad_linear(t,y,A_num,B_num,K,ref),tspan,init_cond);

    % Landing
    Z_mat = states(:,12)*ones(1,12); % make matrix the size of states but each column is z
    States = reshape(states(Z_mat >= 0),[],12);
    Z_mat_lin = states_lin(:,12)*ones(1,12); % make matrix the size of states but each column is z
    States_lin = reshape(states_lin(Z_mat_lin >= 0),[],12);

    x = States(:,10);   x_lin = States_lin(:,10);
    y = States(:,11);   y_lin = States_lin(:,11);
    z = States(:,12);   z_lin = States_lin(:,12);
    phi = States(:,1);  phi_lin = States_lin(:,1);
    th = States(:,2);   th_lin = States_lin(:,2);
    psi = States(:,3);  psi_lin = States_lin(:,3);

    % 2D Visualization
    subplot(rows,2,2*i-1)
    plot_traj(x_lin,z_lin,th_lin,30,LABEL{i,1},LABEL{i,2},LABEL{i,3});
    subplot(rows,2,2*i)
    plot_traj(x,z,th,30,LABEL{i,4},LABEL{i,5},LABEL{i,6});
end
