% Wozniak, Michael P
% AAE568 Final Project
clc; clear; close all;

% Params
n = 0.00113;%mean motion

% this combination works for docking but keeping looks off
Ts = 1;% sampling time (lower as much as possible by compute constraints)
T_sim = 1200; % tune this

% this combination works but is an unreasonable sampling time
Ts = 5;% sampling time (lower as much as possible by compute constraints)
T_sim = 120; % tune this

N = 12;% mpc window
u_limit = 0.05;% .05ms^-2, very conserivative, could be implemented by cold gas thrusters
v_max_dock = 0.01;% based on historical reserach
d_min = 8;% collision avoidance parameter

[Ad, Bd] = CW_dynamics(n, Ts); % fetch discretized state equations based on sampling time and orbit param(s)
K = -dlqr(Ad, Bd, eye(6), eye(3));  % sign flip -dlqr is necessary based on stability of A+BK as set forth in
% page 4, Mammarella
Ak = Ad + Bd*K;% page 4, Mammarella

% Use Tube for tightening constraints
Wbound = 0.0012; % tihs is very tunable but must be really small to get good results
Sinf_x = 0; %  robust positively invariant set Sk for state
Sinf_u = 0; %  robust positively invariant set Sk for control
%approx first 50 (can tune this but 50 is returning robust results)
for j = 0:50
    Sinf_x = Sinf_x + norm(Ak^j, inf)* Wbound;
    Sinf_u = Sinf_u + norm(K*Ak^j, inf)*Wbound;
end
v_lim = u_limit-Sinf_u;
z_lim = 500-Sinf_x; % levying this to keep within a box of 500m just to keep satellites in vicinity of target

% Initial Conditions
d = 250;
x_actual = [d, d, d, 0.2, 0.2, 0.2;
    d, d, -d, -0.2, 0.2, -0.2;
    d,-d, d, 0.2, 0.2, 0.2;
    d,-d, -d, -0.2, 0.2, -0.2;
    -d, d, d, 0.2, 0.2, 0.2;
    -d, d, -d, -0.2, 0.2, -0.2;
    -d,-d, d, 0.2, 0.2, 0.2;
    -d,-d, -d, -0.2, 0.2, -0.2]';
z_nom = x_actual;
x_initial = x_actual;%setpoint=initial point
x_hist = zeros(6, 8, T_sim);%state histry
u_hist = zeros(3, 8, T_sim); % control history
opt = mpcActiveSetOptions;

iA_dock = []; % will be logical on first step
iA_sk = cell(1,6); % logical flags for each observer agent

fprintf('Simulating TRMPC1 & TRMPC2:\n');
activeConsHist = [];
bestActive = 10000;

qGain = [1; 2; 5; 10; 25]';
rGain = [1; 2; 5; 10; 25]';

for qGain = [2.15; 2.20; 2.25; 2.30; 2.35]'
    for rGain = [1.65; 1.7; 1.75; 1.8; 1.85]'
        % Reset states and solver flags for every new gain pair
        z_nom = x_initial;
        x_actual = x_initial;
        iA_dock = [];
        numActive = 0;
        for t = 1:T_sim
            % TRMPC1
            A_dock = blkdiag(Ad, Ad);
            B_dock = blkdiag(Bd, Bd);
            p_rel = z_nom(1:3, 1) - z_nom(1:3, 2);
            [L_d, f_d, E_d, W_d, Aeq_d, beq_d] = prepDockingQP(...
                A_dock, B_dock, N, qGain, rGain, v_lim, z_lim, p_rel, d_min, v_max_dock, z_nom(:,1:2));
            % check if iA_dock is logical
            if t == 1
                iA_dock = false(size(W_d));
            end
            [vseq_d, status_d, iA_dock] = mpcActiveSetSolver(L_d, f_d, E_d, W_d, Aeq_d, beq_d, iA_dock, opt);
            if status_d <= 0
                vseq_d = zeros(72,1);% warning: temporarily storing magic number 72, may break if changing window size
            end
            v_dock = reshape(vseq_d(1:6), 3, 2);

            % Propagate system (with consideration of constraints)
            V_all = v_dock;
            for i = 1:2
                u_actual = V_all(:,i) + K*(x_actual(:,i) - z_nom(:,i));% see paper, true control
                u_actual = max(min(u_actual, u_limit), -u_limit);%hold to u constraint
                u_actual = max(min(u_actual, v_lim), -v_lim);% tube-form, hold to v lim
                if any(abs(abs(u_actual) - v_lim) < 1e-6)
                    numActive = numActive + 1;
                end

                u_hist(:, i, t) = u_actual; % save control
                z_nom(:,i) = Ad*z_nom(:,i) + Bd*V_all(:,i);%update nominal state
                x_actual(:,i) = Ad*x_actual(:,i) + Bd*u_actual + (rand(6,1)-0.5)*Wbound;% propagate actual (with noise)
                x_hist(:, i, t) = x_actual(:,i);% store history
            end
        end % End of T_sim loop

        activeConsHist = [activeConsHist; numActive];
        if numActive < bestActive
            bestActive = numActive;
            qOpt = qGain;
            rOpt = rGain;
        end
        fprintf("Q=%d, R=%d, Total Active Steps=%d\n", qGain, rGain, numActive)
    end
end

fprintf("Min Timesteps of Active Cosntraints: %d\n", bestActive)
fprintf("Q*: %d\n", qOpt)
fprintf("R*: %d\n", rOpt)

function [L, f, E, W, Aeq, beq] = prepDockingQP(A, B, N, qGain, rGain, v_lim, z_lim, p_rel, d_min, v_max, z0)
nx = 12; % 12 elements in state (6 state vec x 2 agents)
nu = 6; % 6 elements in control (ux, uy, uz (or u1/2/3 in paper) and 2 agents)
[H, G] = getHG(A, B, N);
Q = eye(nx*N)*qGain; % tune for weight on state
R = eye(nu*N)*rGain;% tune for weight on control
L = 2*(G'*Q*G + R);
f = 2*(G'*Q*H) * z0(:);

Ev = [eye(nu*N); -eye(nu*N)]; % cpnstraint vec (Across window)
Wv = ones(2*nu*N, 1) * v_lim; % constraint vec (Across window)

Cv = [zeros(3,3), eye(3), zeros(3,6); zeros(3,9), eye(3)]; % extract vel
GN = G(end-nx+1:end, :);
HN = H(end-nx+1:end, :);%,HG across window
Ev_max = [Cv*GN; -Cv*GN];
Wv_max = [ones(6,1)*v_max - Cv*HN*z0(:);
    ones(6,1)*v_max + Cv*HN*z0(:)]; % constraint vec

G_pos = G(1:3, :) - G(7:9, :);
H_pos = H(1:3, :) - H(7:9, :);
E_coll = -p_rel' * G_pos;
W_coll = -d_min*norm(p_rel) + p_rel' * H_pos * z0(:);

E = [Ev; Ev_max; E_coll];
W = [Wv; Wv_max; W_coll];
Aeq = zeros(0, nu*N);
beq = zeros(0, 1);% active set solver particular about Aeq, beq
end
function [L, F, E, W, Aeq, beq] = prepSKQP(A, B, N, v_lim, z_lim)
nx = 6;
nu = 3;
[H, G] = getHG(A, B, N);
Q = eye(nx*N);
R = eye(nu*N)*10;
L = 2*(G'*Q*G + R);
F = 2*(G'*Q*H);
E = [eye(nu*N); -eye(nu*N)];
W = ones(2*nu*N, 1) * v_lim;
Aeq = zeros(0, nu*N);
beq = zeros(0, 1);
end
function [H, G] = getHG(A, B, N)
nx = size(A,1);
nu = size(B,2);
H = [];
G = [];
for i = 1:N
    H = [H; A^i];
    row = [];
    for j = 1:N
        if j <= i
            row = [row, A^(i-j)*B];
        else
            row = [row, zeros(nx,nu)];
        end
    end
    G = [G; row];
end
end
function [Ad, Bd] = CW_dynamics(n, Ts)
% continuous versions
Ac = [zeros(3), eye(3);
    3*n^2, 0, 0, 0, 2*n, 0; 0, 0, 0, -2*n, 0, 0;
    0, 0, -n^2, 0, 0, 0];
Bc = [zeros(3);
    eye(3)];
sys = c2d(ss(Ac, Bc, eye(6), 0), Ts);% discretize
% discrete vresions
Ad = sys.A;
Bd = sys.B;
end