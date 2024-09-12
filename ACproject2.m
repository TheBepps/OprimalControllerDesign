clear
close all
clc

%% Data
g = 9.81;               %Gravitational constant[m/s^2]
md = 0.5;               %Drone mass[kg]
Iz = 0.05;              %Drone inertia[kg m^2]

%Control inputs
f = sdpvar(1,1);        %Thrusting force  
tauc= sdpvar(1,1);      %Control torque

%Error coordinates 
ez = sdpvar(1,1);
epsi = sdpvar(1,1);

%Dynamics
vz = sdpvar(1,1);       %z velocity
omegaz = sdpvar(1,1);   %z angular speed

% To remove the −g term in the linear acceleration dynamics you can define an auxiliary 
% input ˆf = f−mdg and write the State Space representation using as input vector u = [ˆf τ]^T ∈ R2
fhat = f-md*g;

%% State space rappresentation
x = [ez;vz;epsi;omegaz];  %Dynamics matrix
u = [fhat;tauc];          %Input matrix
z = [ez;epsi];            %Output matrix

A = [ 0   1   0   0;
    0   0   0   0;
    0   0   0   1;
    0   0   0   0];
B = [ 0       0;
    1/md    0;
    0       0;
    0       1/Iz];
C = [ 1   0   0   0;
    0   0   1   0];
D = [ 0   0;
    0   0];

[n,~] = size(A);
[~,p] = size(B);

%%  1)Is the system controllable? 
rank_R = rank(ctrb(A,B));
if(rank_R == n)
    disp('Controllability matrix is full rank -> system is controllable')
else
    disp('System is not controllable')
end

%%  2)Design of a full-state static feedback controller 
yalmip('clear');
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

W=sdpvar(n,n);
X=sdpvar(p,n);
fake_zero=1e-3;
alpha=1;
k=sdpvar(1,1);

constr=[k==200;
        W>=fake_zero*eye(n);
        A*W+B*X+(A*W+B*X)'<=-2*alpha*W;
        [k*eye(n), X';
         X, k*eye(p)] >= fake_zero*eye(n+p)];

solvesdp(constr, k, opts);

W_a=value(W);
X_a=value(X);
P_a=inv(W_a);
K_a=X_a*inv(W_a);

[primal_res, dual_res] = check(constr);

if( all([primal_res; dual_res] > -fake_zero) )
    disp('-----------------------------------------------------------------');
    disp('Controller imposing convergence rate and minimizing actuator effort');
    disp('Problem is feasible -> linearized system is asymptotically stable');
    disp('Gain matrix');
    disp(K_a);
    disp('Real part of closed-loop eigenvalues');
    disp(mat2str(real(eig(A+B*K_a)')));
else
    disp('Problem is infeasible');
end

%% 3)State space representation with perturbations
AA = [  0   1   0   0   0   0;
        0   0   1   0   0   0;
        0   0   0   0   0   0;
        0   0   0   0   1   0;
        0   0   0   0   0   1;
        0   0   0   0   0   0];

BB = [  0     0;
        0     0;
        1/md  0;
        0     0;
        0     0;
        0     1/Iz];

CC = [  0   1   0   0   0   0;
        0   0   0   0   1   0];

DD = [  0   0;
        0   0];

EE = [  0    0;
        0    0;
        1/md 0;
        0    0;
        0    0;
        0    1/Iz];

FF = [  0   0;
        0   0];

[n_p,~] = size(AA); %[n x n]
[~,p_p] = size(BB); %[n x p]
[m,~] = size(CC); %[m x n]
[~,d] = size(FF); %[m x d]

%% Optimality curve (3.1)

yalmip('clear')
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

W_p = sdpvar(n_p,n_p);
X_p = sdpvar(p_p,n_p);
rho = sdpvar(1,1);
gamma = sdpvar(1,1);
fake_zero = 1e-3;

gamma_opt = [];
k_opt = [];

% -- k in the interval [10, 1400]
for k = 10:10:1400
    constr = [W_p >= rho*eye(n_p);

              rho >= fake_zero;

              [k*rho*eye(n_p),        X_p';
                  X_p,          k*rho*eye(p_p)] >= fake_zero*eye(n_p+p_p);

              [AA*W_p + BB*X_p + (AA*W_p + BB*X_p)',      EE,     (CC*W_p + DD*X_p)';
                                EE',                   -gamma*eye(d),         FF';
                      CC*W_p + DD*X_p,                      FF,        -gamma*eye(m)] <= -fake_zero*eye(n_p+d+m)];
                  
                  
     solvesdp(constr, gamma, opts);
     
     [primal_res, dual_res] = check(constr);
     
     if( all([primal_res; dual_res] > -fake_zero) )
         gamma_opt(end+1) = value(gamma);
         k_opt(end+1) = k;
     else
         disp('-----------------------------------------------------------------');
         disp('Problem is infeasible'); 
     end
    
end

% -- Plot optimality curve
figure(1)
area(k_opt, gamma_opt, 'FaceColor',[0.3010 0.7450 0.9330], 'EdgeColor',[0 0.4470 0.7410]);
xlabel('Bound on |K|');
ylabel('Minimum value of \gamma');
axis([0 k_opt(end) 0 gamma_opt(1)+1/4]);

%% gamma < 0.1, alpha = 1, |K| < 200 (3.2)
% -- Impose alpha = 1, gamma < 0.1 (for example gamma = 0.09) 
alpha = 1;
k_des = 200;

constr = [  W_p >= rho*eye(n_p);

            rho >= fake_zero;

          [k_des*rho*eye(n_p),         X_p';
                  X_p,          k_des*rho*eye(p_p)] >= fake_zero*eye(n_p+p_p);

          [AA*W_p + BB*X_p + (AA*W_p + BB*X_p)',      EE,     (CC*W_p + DD*X_p)';
                            EE',                   -gamma*eye(d),         FF';
                  CC*W_p + DD*X_p,                      FF,        -gamma*eye(m)] <= -fake_zero*eye(n_p+d+m);

           AA*W_p + BB*X_p + (AA*W_p + BB*X_p)' <= -2*alpha*W_p;

           gamma <= 0.1 - fake_zero];

% -- Solve the problem
sol = solvesdp(constr,[],opts);

[primal_res, dual_res] = check(constr);
disp('-----------------------------------------------------------------');
if( all([primal_res; dual_res] > -fake_zero) )
    disp('The problem with bounded parameters is feasible.')
else
    disp('The problem with bounded parameters is infeasible.');
end

%% Colleague solution with k = 1200 (3.3)
k_col = 1200;

constr = [W_p >= rho*eye(n_p);

          rho >= fake_zero;

          [k_col*rho*eye(n_p),         X_p';
                  X_p,          k_col*rho*eye(p_p)] >= fake_zero*eye(n_p+p_p);

          [AA*W_p + BB*X_p + (AA*W_p + BB*X_p)',      EE,     (CC*W_p + DD*X_p)';
                            EE',                   -gamma*eye(d),         FF';
                  CC*W_p + DD*X_p,                      FF,        -gamma*eye(m)] <= -fake_zero*eye(n_p+d+m);

           AA*W_p + BB*X_p + (AA*W_p + BB*X_p)' <= -2*alpha*W_p;

           gamma <= 0.1 - fake_zero];

% -- Solve the problem
solvesdp(constr, gamma, opts);

[primal_res, dual_res] = check(constr);

if( all([primal_res; dual_res] > -fake_zero) )
    disp('The problem with colligue parameters is feasible.')
else
    disp('The problem with colligue parameters is infeasible.');
end

% Extract values
W_val = value(W_p);
X_val = value(X_p);
K_val = X_val / W_val;
gamma = value(gamma);
K_norm = norm(K_val);


disp('The gain matrix K is:')
disp(value(K_val));
disp('The norm of K is:')
disp(value(K_norm));
disp('The real part of the closed loop eigenvalues is:')
disp(value(real(eig(AA+BB*K_val)')));
disp('gamma is:')
disp(value(gamma));

%% Simulation of the system from initial conditios (3.4)
t_f = 15;
step_size = 1e-3;
t_sim = 0:step_size:t_f;

x_0 = [0, 1, 0.1, 0, 0.5 0]';

w = [2, 0.2]';

[t_ode_CL, x_ode_CL] = ode45( @(t,x)(AA*x + BB*K_val*x + EE*w), t_sim, x_0 ); 

z_CL = CC*x_ode_CL';

figure(3);
plot(t_sim, z_CL(1,:),'color','b','linewidth',1)
hold on
plot(t_sim, z_CL(2,:),'color','g','linewidth',1)
hold off
title('Feedback-control z and psi errors')
xlabel('Time [s]');
ylabel('Errors');
legend('z-error','psi-error');