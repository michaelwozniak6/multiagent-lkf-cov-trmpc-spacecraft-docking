%% Wozniak, Michael P
% AAE568 Final Project 
clc; clear; close all;

%% First - Optimal Estimation for IOD over 1200s
fprintf('Simulating Multi-Agent Linear KF\n');
num_agents = 8;
n = 0.00113;
t0 = 0;
tf = 1200; % consistent with COV implementation, reasonable time-horizon for this mission
tvec = linspace(t0, tf, tf);
M = [eye(3), zeros(3)]; % because position only
W_inv = eye(3) * (2.25^2); % measurement noise
all_norm_pos_errors = cell(num_agents, 1); % will append respective errors for each agent
all_final_pos = cell(num_agents, 1); % will append final positions for each agent
% Define reference trajectory (init guess)
d = 250; % distance parameter for initial spacing
x0_ref = -.5*[d, d,  d, 0.2, 0.2, 0.2]'; % initial guess: in the vicinity of C_1
x0actual = -[d, d,  d, 0.2, 0.2, 0.2]'; % agent 1 is the coordinate frame oriign for KF purposes, so target is actually initially at the negated x0_symmetric for C_1
x0_symmetric =[d, d,  d, 0.2, 0.2, 0.2;
     d, d, -d, -0.2, 0.2, -0.2;
     d,-d,  d, 0.2, 0.2, 0.2;
     d,-d, -d, -0.2, 0.2, -0.2;
    -d, d,  d, 0.2, 0.2, 0.2;
    -d, d, -d, -0.2, 0.2, -0.2;
    -d,-d,  d, 0.2, 0.2, 0.2;
    -d,-d, -d, -0.2, 0.2, -0.2]'; % symmetric pos, different initial velocities

% necessary for the dynamics (and translating to C_1 LVLH in KF)
Re = 6378137;%earth radius m
h = 500000;%500km altitude
Rorb = Re + h;

for agentindex = 1:num_agents
    % chaser spacecraft formation
    offsets = x0_symmetric - x0_symmetric(:,1); %shift all init pos so agent 1 is at the origin (will subsequently use DCM to translate measuremnets)
    % Simulate measurements
    yM = cell(agentindex, 1);
    for k = 1:agentindex
        for i = 1:length(tvec)
            Phi_true = cwSTM(tvec(i), n); % Mapping t0 to ti
            x_target_ti = Phi_true * x0actual;
            % Measurement is target at ti minus agent offset
            dy = offsets(2,k);
            theta=dy/Rorb;
            RC1Ck=[cos(theta), sin(theta), 0;-sin(theta),cos(theta),0;0, 0, 1];%DCM for translating measurements to C1
            initialOffsetTrue = x0_symmetric(:,k) - x0_symmetric(:,1);%start update for propagating offset each step
            currentOffsetTrue = Phi_true*initialOffsetTrue;
            dyTrue = currentOffsetTrue(2);
            thetaTrue = dyTrue/Rorb;
            RC1CkTrue = [cos(thetaTrue), sin(thetaTrue), 0; -sin(thetaTrue), cos(thetaTrue), 0; 0, 0, 1];
            yM{k}(:,i) = RC1CkTrue * (M * x_target_ti - currentOffsetTrue(1:3)) + 2.25*randn(3,1);% end update propagaating each step
        end
    end

    % Init LKF
    x_plus = x0_ref; % x_{i-1}^+
    P_plus = eye(6) * 1e8; % initial covariance  (high init guess)
    xhist = zeros(6, length(tvec));
    for i = 1:length(tvec)
        % STM Mapping - Phi_i is the STM from t_{i-1} to t_i for the time update
        % (for the measurement G, will need Phi from t_0 to t_i)
        dt = 1; % unit step for simplicity
        Phi_step = cwSTM(dt, n);
        Phi_from_t0 = cwSTM(tvec(i), n);
        
        %Reference state at current time
        xi_ref = Phi_from_t0 * x0_ref;
        x_prev_ref = cwSTM(tvec(max(1,i-1)), n) * x0_ref;

        %TIME update (see p352-354 lec notes)
        % x_i^- = x_i^ref + Phi_i(x_{i-1}^+ - x_{i-1}^ref)
        x_minus = xi_ref + Phi_step * (x_plus - x_prev_ref);
        Q = eye(6) * 1e-9; % model uncertainty
        P_minus = Phi_step * P_plus * Phi_step' + Q;
        % MEASUREMENT Update
        for k = 1:agentindex
            initialOffsetK = x0_symmetric(:,k) - x0_symmetric(:,1);
            currentOffsetK = Phi_from_t0*initialOffsetK;
            dy = offsets(2,k);
            theta = dy/Rorb;
            RC1Ck = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
            % Gi = partial g / partial x at reference
            Gi = RC1Ck*M; % DCM must be used for the dynamics
            % z_i: actual measurement
            zi = yM{k}(:,i);
            % g_i(xi_ref): what the measurement would be if state = reference
            % Account for the agent offset in the measurement model
            gi_ref = RC1Ck*(M * xi_ref) -(RC1Ck* currentOffsetK(1:3));% DCM must be used for the dynamics
            %K_i = P_i^- Gi' (W_inv + Gi P_i^- Gi')^-1
            Si = Gi * P_minus * Gi' + W_inv;
            Ki = P_minus * Gi' / Si;
            %x_i^+ = x_i^- + Ki(zi - gi_ref - Gi(x_minus - xi_ref))
            innov = zi - gi_ref - Gi * (x_minus - xi_ref);
            x_plus = x_minus + Ki * innov;
            % P_i^+ = (I - Ki Gi) * P_i^-
            P_plus = (eye(6) - Ki * Gi) * P_minus;
            % Plus for this timestep is minus for next time step
            x_minus = x_plus;
            P_minus = P_plus;
        end
        xhist(:,i) = x_plus;% Store the estimate of the state at time t_i
    end
    % Normalized Position Error to quantify success (filter tracks state at
    % time t so compare to xtrue at time t)
    err_at_t = zeros(1, length(tvec));
    for i = 1:length(tvec)
        xtrue_i = cwSTM(tvec(i), n) * x0actual;
        err_at_t(i) = norm(xhist(1:3,i) - xtrue_i(1:3)) / (norm(xtrue_i(1:3)) + 1e-3);%compute with conventional error form
    end
    all_norm_pos_errors{agentindex} = err_at_t;
    all_final_pos{agentindex} = xhist(1:6,end);
end

% Fig: 3 subplots of norm of position error estimate with time for
% different # agnets
figure;
kprime=1;
for k = [1, 4, 8]
    subplot(1,3,kprime);
    plot(tvec, all_norm_pos_errors{k}, 'LineWidth', 1.5);
    grid on; 
    title(sprintf('%d Agents', k));
    ylabel('Norm of Position Error (m)')
    xlabel('Time (s)');
    kprime=kprime+1;
end
sgtitle('Error in Estimation for Varying Agents')
% Fig; shows benefit of increasing # agents
means = [];
figure;
for it= 1:num_agents
    means=[means; mean(all_norm_pos_errors{it})];
end
plot(1:num_agents, means,'o--','LineWidth',2, 'MarkerFaceColor','blue', 'Color','blue');
grid on
xlabel('Number of Agents')
ylabel('Mean of Norm Position Error (m)')
title('Mean of Norm Pos. Error with Num. Agents')

x0_symmetric=x0_symmetric';
% Fig; shows problem formulation via initial layout
figure;
hold on; 
grid on;
axis equal;
% Target at origin (this is in x0_symmetric back before putting LVLH on C1
% and estimating wrt that - just for illustraiton)
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'Target Spacecraft');
% 8 Chaser Agents
scatter3(x0_symmetric(:,1), x0_symmetric(:,2), x0_symmetric(:,3), 80, 'filled','MarkerFaceColor', 'b',  'DisplayName', 'Chaser Spacecraft');
% Line of Sight to help visualize chaser positions
for i = 1:size(x0_symmetric, 1)
    plot3([x0_symmetric(i,1), 0], [x0_symmetric(i,2), 0], [x0_symmetric(i,3), 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    text(x0_symmetric(i,1)+5, x0_symmetric(i,2)+5, x0_symmetric(i,3)+5, sprintf('C_%d', i), 'FontSize', 10, 'FontWeight', 'bold');
end
% Coordinate Frame Vectors at Agent 1 C1
% Position of Agent 1
C1 = x0_symmetric(1, 1:3); 
vScale = 100; % Length of the coordinate axis arrows for visibility
q1 = quiver3(C1(1), C1(2), C1(3), vScale, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'x (Radial)');
q2 = quiver3(C1(1), C1(2), C1(3), 0, vScale, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'y (In-track)');
q3 = quiver3(C1(1), C1(2), C1(3), 0, 0, vScale, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'z (Cross-track)');
% 5. Formatting
xlabel('x (m)'); 
ylabel('y (m)'); 
zlabel('z (m)');
title('Initial Relative Positions and C1 LVLH Coordinate Frame');
legend('Location', 'northeastoutside');
view(3);
% Geometry of C1 Estimate
target_pos = [0, 0, 0];
C1_pos = x0_symmetric(1, 1:3);
% direction vector y1
y1_vec = target_pos - C1_pos; 
q_obs = quiver3(C1_pos(1), C1_pos(2), C1_pos(3), y1_vec(1), y1_vec(2), y1_vec(3), 'AutoScale', 'off', 'Color', 'cyan', 'LineWidth', 2, 'DisplayName', 'Observation y_1');
midpoint = C1_pos + 0.5 * y1_vec;
text(midpoint(1), midpoint(2)+15, midpoint(3), '$y_1(t_i)$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

%% Second - Online Calculus of Variations to Determine Agents which would be
% most optimal to dock with target for minimum cost (min control effort)
% x0 of task selection is xfinal of LKF (assuming it took 1200s for debris
% cataloging)
fprintf('Optimally Selecting Agents via Calculus of Variations\n');
xTargetRelC1_tf = all_final_pos{8};  % case of 8 observers is most accurate
x0TaskSelection = zeros(num_agents, 6); % will have to reconstruct our init state matrix
% (after performing transformaitons to get all in same frame)
Phi_tf = cwSTM(tf, n); %Propagate 1200s
for i = 1:num_agents
    % Calculate agnet i's current offset relative to C1 at tf considering initial table values from x0_symmetric
    initial_offset_rel_C1 = x0_symmetric(i,:)' - x0_symmetric(1,:)';
    current_offset_rel_C1 = Phi_tf * initial_offset_rel_C1;
    % Agent i relative to Target = (Agent i rel C1) - (Target rel C1)
    % place Target at the Origin (0,0,0) for the next phase
    x0TaskSelection(i, :) = (current_offset_rel_C1 - xTargetRelC1_tf)';
end

Jvals = [];
for agent = 1:numel(x0TaskSelection(:,1))
    J = task_assignment(x0TaskSelection(agent,:)', agent);
    Jvals = [Jvals; J];
    fprintf("Agent %d: J*=%.8f\n", agent, J)
end
[minCosts, chosenAgents]=mink(Jvals,2);
fprintf("Selected Agents:\n")
fprintf("Agent %d: J*=%.8f\n", chosenAgents(1), minCosts(1))
fprintf("Agent %d: J*=%.8f\n", chosenAgents(2), minCosts(2))

% Plot initial positions (after surveying period)
figure;
hold on; 
grid on;
axis equal;
% Target at origin
plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'Target');
% Plot 8 Chaser Agents
scatter3(x0TaskSelection(:,1), x0TaskSelection(:,2), x0TaskSelection(:,3), 80, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Chaser Agents');
for i = 1:size(x0TaskSelection, 1)
    text(x0TaskSelection(i,1)+5, x0TaskSelection(i,2)+5, x0TaskSelection(i,3)+5, sprintf('C_%d', i), 'FontSize', 10, 'FontWeight', 'bold');
end
for i = 1:size(x0_symmetric, 1)
    plot3([x0TaskSelection(i,1), 0], [x0TaskSelection(i,2), 0], [x0TaskSelection(i,3), 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    text(x0TaskSelection(i,1)+5, x0TaskSelection(i,2)+5, x0TaskSelection(i,3)+5, sprintf('C_%d', i), 'FontSize', 10, 'FontWeight', 'bold');
end
xlabel('x (m)'); 
ylabel('y (m)'); 
zlabel('z (m)');
title('Relative Positions of Chaser Spacecraft after Debris Cataloging Period');
legend('Location', 'northeastoutside');
view(3);

%% Tube-MPC
% Params
% implemented TMPC based on 1,2 agents being optimal; now 7,5 are most
% optimal. have to do a mapping / Change of vars
allIndices = 1:8;
observers = setdiff(allIndices, chosenAgents);
newOrder = [chosenAgents', observers];

n = 0.00113;%mean motion

Ts = 5;% sampling time (lower as much as possible by compute constraints)
T_sim = 120*2; % tune this

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
z_lim = 500-Sinf_x; % this would keep within a box of 500m 

% Initial Conditions
d = 250; 
x_actual = x0TaskSelection(newOrder,:)';
z_nom = x_actual; 
x_initial = x0_symmetric(newOrder,:)';%setpoint=initial point (not just from task selection but from entire time-history)
x_hist = zeros(6, 8, T_sim);%state histry
u_hist = zeros(3, 8, T_sim); % control history
opt = mpcActiveSetOptions;

iA_dock = []; % will be logical on first step
iA_sk = cell(1,6); % logical flags for each observer agent

fprintf('Simulating TRMPC1 & TRMPC2\n');
for t = 1:T_sim
    % TRMPC1
    A_dock = blkdiag(Ad, Ad); 
    B_dock = blkdiag(Bd, Bd);
    p_rel = z_nom(1:3, 1) - z_nom(1:3, 2); 
    [L_d, f_d, E_d, W_d, Aeq_d, beq_d] = prepDockingQP(...
        A_dock, B_dock, N, v_lim, z_lim, p_rel, d_min, v_max_dock, z_nom(:,1:2));
    % initialize iA to logical false
    if t == 1
        iA_dock = false(size(W_d));
    end
    [vseq_d, status_d, iA_dock] = mpcActiveSetSolver(L_d, f_d, E_d, W_d, Aeq_d, beq_d, iA_dock, opt);
    if status_d <= 0
        vseq_d = zeros(72,1);% warning: temporarily storing magic number 72, may break if changing window size, just state * window (12*6)
    end
    v_dock = reshape(vseq_d(1:6), 3, 2);

    % TRMPC2
    [L_s, F_s, E_s, W_s, Aeq_s, beq_s, H_s] = prepSKQP(Ad, Bd, N, v_lim, z_lim);
    v_sk = zeros(3, 6);
    for i = 1:6
        z_k = z_nom(:,i+2);
        x_ref = x_initial(:,i+2);
        
        % build stacked reference over horizon
        X_ref = repmat(x_ref, N, 1);
        
        f_s = F_s * (H_s*z_k - X_ref);
        % initialize iA to logical false
        if t == 1
            iA_sk{i} = false(size(W_s));
        end
        [v_seq_s, status_s, iA_sk{i}] = mpcActiveSetSolver(L_s, f_s, E_s, W_s, Aeq_s, beq_s, iA_sk{i}, opt);
        if status_s <= 0
            v_seq_s = zeros(36,1);% warning: storing magic number, may break if you change window size, just window size * num contorl (12*3)
        end
        v_sk(:, i) = v_seq_s(1:3);
    end

    % Propagate system (with consideration of constraints)
    V_all = [v_dock, v_sk];
    x_est = kfStep(x_actual, x0TaskSelection);
    for i = 1:8
        u_actual = V_all(:,i) + K*(x_actual(:,i) - z_nom(:,i));% see mammarella paper
        u_actual = max(min(u_actual, u_limit), -u_limit);% then efficienctly hold to min or max where needed
        u_hist(:, i, t) = u_actual; % save control
        z_nom(:,i) = Ad*z_nom(:,i) + Bd*V_all(:,i);%update nominal state
        % x_actual(:,i) = Ad*x_actual(:,i) + Bd*u_actual + (rand(6,1)-0.5)*Wbound;% propagate actual (with noise)
        g = nonlinear_perturbation(x_actual(:,i),n);
        x_actual(:,i) = Ad*x_actual(:,i) + Bd*u_actual + Bd*g + (rand(6,1)-0.5);% propagate actual (with noise)
        x_hist(:, i, t) = x_actual(:,i);% store history
    end
end

% Time-History and Tube Visualization
figure;
hold on; 
grid on; 
axis equal; 
view(3);
colors = lines(8);
[sx, sy, sz] = sphere(20);
sx = sx * Sinf_x; 
sy = sy * Sinf_x; 
sz = sz * Sinf_x;

% Docking agents
for i = 1:2
    actual_pos = squeeze(x_hist(1:3, i, :));
    
    % Initial point (marker necessary so it's visible)
    plot3(actual_pos(1,1), actual_pos(2,1), actual_pos(3,1), '.', 'Color', colors(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off');
    
    % Actual Traj
    % plot3(actual_pos(1,:), actual_pos(2,:), actual_pos(3,:), 'Color', colors(i,:), ...
    %     'LineWidth', 2, 'DisplayName', sprintf('Actual State xk for C%d', i));
    st = '-'; if i > 2, st = '--'; end
    plot3(actual_pos(1,:), actual_pos(2,:), actual_pos(3,:), 'LineStyle', st, ...
        'Color', [colors(i,:), 0.4], 'LineWidth', 2, ... 
        'DisplayName', sprintf('Actual State xk for C%d', newOrder(i)));
    
    % Nominal Spine
    % plot3(actual_pos(1,:), actual_pos(2,:), actual_pos(3,:), ':', 'Color', colors(i,:)*0.5, ...
    %     'LineWidth', 1, 'DisplayName', sprintf('Nominal Spine zk for C%d', i));
    plot3(actual_pos(1,:), actual_pos(2,:), actual_pos(3,:), ':', ...
        'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Nominal Spine zk for C%d', newOrder(i)));
end

% Observer Agents
for ch = 3:8
    actual_pos = squeeze(x_hist(1:3, ch, :));
    % Initial position marker
    plot3(x_initial(1, ch), x_initial(2, ch), x_initial(3, ch), '.', 'Color', colors(ch,:), ...
        'MarkerSize', 25, 'DisplayName', sprintf('Station C%d (Hold Pos)', newOrder(ch)));
    % Trajectory (Dashed, if visible)
    plot3(actual_pos(1,:), actual_pos(2,:), actual_pos(3,:), '--', 'Color', colors(ch,:), ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end
% Formatting
plot3(0,0,0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'Target');
xlabel('x (m)'); 
ylabel('y (m)'); 
zlabel('z (m)');
title('Full Time-History of Chasers C1-C8');
% debug strat to get all observers visible
axis tight; 
lims = [xlim; ylim; zlim];
xlim([lims(1,1)-20, lims(1,2)+20]);
ylim([lims(2,1)-20, lims(2,2)+20]);
zlim([lims(3,1)-20, lims(3,2)+20]);
% Legend

legend('Interpreter', 'latex', 'Location', 'northeastoutside');

% Docking Agents history
time = (0:T_sim-1) * Ts;
figure('Color', 'w');
subplot(2,1,1); 
hold on; 
grid on;
for i = 1:2
    plot(time, sqrt(sum(squeeze(x_hist(1:3, i, :)).^2, 1)), 'LineWidth', 1.5);
end
ylabel('Pos Norm [m]'); 
title('Docking: Relative Position to Target');
% legend('C7', 'C5')
legend(sprintf('C%d', newOrder(1)), sprintf('C%d', newOrder(2)))
subplot(2,1,2); 
hold on; 
grid on;
for i = 1:2
    plot(time, sqrt(sum(squeeze(u_hist(:, i, :)).^2, 1)), 'LineWidth', 1.5);
end
ylabel('Control [m/s^2]'); 
xlabel('Time [s]'); 
legend(sprintf('C%d', newOrder(1)), sprintf('C%d', newOrder(2)))
title('Docking: Control Time History');

% Station Keeping history
figure('Color', 'w');
subplot(2,1,1); 
hold on; 
grid on;
for i = 3:8
    plot(time, sqrt(sum(squeeze(x_hist(1:3, i, :)).^2, 1)), 'LineWidth', 1.5);
end
ylabel('Pos Norm [m]'); 
% legend('C1', 'C2', 'C3', 'C4', 'C6', 'C8');
sk_labels = arrayfun(@(id) sprintf('C%d', id), newOrder(3:8), 'UniformOutput', false);
legend(sk_labels);
title('Keeping: Position Norm from Target');
subplot(2,1,2); 
hold on; 
grid on;
for i = 3:8
    plot(time, sqrt(sum(squeeze(u_hist(:, i, :)).^2, 1)), 'LineWidth', 1.5);
end
ylabel('Control [m/s^2]'); 
xlabel('Time [s]'); 
title('Keeping: Control Norm Time History');
% legend('C1', 'C2', 'C3', 'C4', 'C6', 'C8');
sk_labels = arrayfun(@(id) sprintf('C%d', id), newOrder(3:8), 'UniformOutput', false);
legend(sk_labels);

%% Function definitions
function Phi = cwSTM(t, n)
    nt = n*t; 
    s = sin(nt); 
    c = cos(nt);
    Phi_rr = [4-3*c, 0, 0; 6*(s-nt), 1, 0; 0, 0, c];
    Phi_rv = [(1/n)*s, (2/n)*(1-c), 0; (2/n)*(c-1), (1/n)*(4*s-3*nt), 0; 0, 0, (1/n)*s];
    Phi_vr = [3*n*s, 0, 0; 6*n*(c-1), 0, 0; 0, 0, -n*s];
    Phi_vv = [c, 2*s, 0; -2*s, 4*c-3, 0; 0, 0, c];
    Phi = [Phi_rr Phi_rv; Phi_vr Phi_vv];
end


function Jstar = task_assignment(x0, agent)
n = 0.0011;
t0 = 0;
tf = 1200;

xf = [0; 0; 0; 0; 0; 0];
% solinit = bvpinit(linspace(t0, tf, 100), zeros(13, 1));
% solinit = bvpinit(linspace(t0, tf, 100), [zeros(6, 1); ones(6,1);1]);
solinit = bvpinit(linspace(t0, tf, 100), [zeros(6, 1); 0.1*ones(6,1);.005]);
options = bvpset('RelTol', 1e-4);
sol = bvp4c(@(t, X) spacecraft_dynamics(t, X, n), @(xa, xb) [xa(1:6)-x0; xb(1:6)-xf; xa(13)-0], solinit, options);
t_hist = sol.x;
x_hist = sol.y(1:6, :);
l_hist = sol.y(7:12,:);
Jstar = sol.y(13,end);
rel_dist = sqrt(x_hist(1,:).^2 + x_hist(2,:).^2 + x_hist(3,:).^2);

umax=0.05; % tune this based on mission-ready requirements
% warning umax defined in two places
u = zeros(numel(t_hist),1);
u_ll = 0;
u_ul = umax;
constraints = zeros(numel(u),1);

for i = 1:length(t_hist)
    mu1=0;
    mu2=0;
    mu3=0;
    u1_unconstrained = (mu1-l_hist(4,i));
    u2_unconstrained = (mu2-l_hist(5,i));
    u3_unconstrained = (mu3-l_hist(6,i));
    u1(i) = u1_unconstrained;
    u2(i) = u2_unconstrained;
    u3(i) = u3_unconstrained;
    % if c_1,c2,c3 are not violated:
    if u1(i) > umax
        mu1 = -l_hist(4,i) - umax;
        u1(i) = umax;
    end
    if u2(i) > umax
        mu2 = -l_hist(5,i) - umax;
        u2(i) = umax;
    end
    if u3(i) > umax
        mu3 = -l_hist(6,i) - umax;
        u3(i) = umax;
    end
end
if (agent == 5) || (agent == 7)
    % hard-coding the plotting of optimal values, will have to update if
    % changing problem form
    figure;
    subplot(2, 1, 1);
    plot(t_hist, rel_dist, 'k', 'LineWidth', 2); hold on;
    grid on;
    xlabel('Time (s)'); ylabel('Norm of Relative Distance (m)');
    title(sprintf('Time-History from COV-Based Task Selection: Optimal Agent %d', agent));

    L_history = sol.y(7:12, :);
    % Derived the optimal control law u = -B' * L
    % Since B' = [zeros(3,3), eye(3)], u is just the negative of the last 3 costates
    u_history = -L_history(4:6, :);
    subplot(2, 1, 2);
    plot(t_hist, u_history, 'LineWidth',1.75)
    grid on;
    xlabel('Time (s)'); ylabel('Optimal Control (m/s^2)');
    title(sprintf('Optimal Control from COV-Based Task Selection: Optimal Agent %d', agent));
    legend('u_1','u_2','u_3')
end
end
function dydt = spacecraft_dynamics_linear(t, y, n)
% n = 0.0011;
x=y(1:6);
L=y(7:12);
lambda1=L(1);
lambda2=L(2);
lambda3=L(3);
lambda4=L(4);
lambda5=L(5);
lambda6=L(6);

umax=0.05; % tune this based on mission-ready requirements
% if c_1,c2,c3 are not violated:
mu1=0;mu2=0;mu3=0; % will overwrite if no Ci violated
u_ll = 0;
u_ul = umax;

u1_unconstrained = (mu1-lambda4);
u2_unconstrained = (mu2-lambda5);
u3_unconstrained = (mu3-lambda6);
u1 = u1_unconstrained;
u2 = u2_unconstrained;
u3 = u3_unconstrained;
% if c_1,c2,c3 are not violated:
if u1 > umax
    mu1 = -lambda4 - umax;
    u1 = umax;
end
if u2 > umax
    mu2 = -lambda5 - umax;
    u2 = umax;
end
if u3 > umax
    mu3 = -lambda6 - umax;
    u3 = umax;
end
u=[u1;u2;u3];
A = [0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    3*n^2, 0, 0, 0, 2*n, 0;
    0, 0, 0, -2*n, 0, 0;
    0, 0, -n^2, 0, 0, 0];
B = [zeros(3); eye(3)];
dx = A * x(1:6) + B * u;
dL1=-3*(n^2)*lambda4;
dL2=0;
dL3=(n^2)*lambda6;
dL4=-(lambda1-2*n*lambda5);
dL5=-(lambda2+2*n*lambda4);
dL6=-lambda3;
dL=[dL1;dL2;dL3;dL4;dL5;dL6];
% dydt = [dx; dL]; former
dydt = [dx; dL; 0.5*u'*u]; % trying to compute cost
end
function dydt = spacecraft_dynamics(t, y, n)
% n = 0.0011;
x=y(1:6);
L=y(7:12);
lambda1=L(1);
lambda2=L(2);
lambda3=L(3);
lambda4=L(4);
lambda5=L(5);
lambda6=L(6);

umax=0.05; % tune this based on mission-ready requirements
% if c_1,c2,c3 are not violated:
mu1=0;
mu2=0;
mu3=0; 
% will overwrite if no Ci violated
u_ll = 0;
u_ul = umax;
u1_unconstrained = (mu1-lambda4);
u2_unconstrained = (mu2-lambda5);
u3_unconstrained = (mu3-lambda6);
u1 = u1_unconstrained;
u2 = u2_unconstrained;
u3 = u3_unconstrained;
% if c_1,c2,c3 are not violated:
if u1 > umax
    mu1 = -lambda4 - umax;
    u1 = umax;
end
if u2 > umax
    mu2 = -lambda5 - umax;
    u2 = umax;
end
if u3 > umax
    mu3 = -lambda6 - umax;
    u3 = umax;
end
u=[u1;u2;u3];

rE = 6378137;
J2 = 1.08263e-3;
mu = 3.986e14;
rTarget = rE + 500000;
A = [0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    3*n^2, 0, 0, 0, 2*n, 0;
    0, 0, 0, -2*n, 0, 0;
    0, 0, -(n^2), 0, 0, 0];

B = [zeros(3); eye(3)];

% nonlinear correction terms
rx = rTarget + x(1);
ry = x(2);
rz = x(3);
R = sqrt(rx^2 + ry^2 + rz^2);
g1 = -(mu/R^3)*rx - 3*n^2*rx;
g2 = -(mu/R^3-n^2)*ry;
g3 = -rz*(mu*rz/R^3-n^2);
r = R;
g = [g1;
     g2;
     g3];

dx = A * x(1:6) + B * u + B * g;
dL1=-3*(n^2)*lambda4;
dL2=0;
dL3=(n^2)*lambda6;
dL4=-(lambda1-2*n*lambda5);
dL5=-(lambda2+2*n*lambda4);
dL6=-lambda3;
dL=[dL1;dL2;dL3;dL4;dL5;dL6];

dydt = [dx; dL; 0.5*u'*u]; % trying to compute cost
% added one more element in dydt that is the integrand
end
function g = nonlinear_perturbation(x, n)
% n = 0.0011;
rE = 6378137;
mu = 3.986e14;
rTarget = rE + 500000;

% nonlinear correction terms
rx = rTarget + x(1);
ry = x(2);
rz = x(3);
R = sqrt(rx^2 + ry^2 + rz^2);
g1 = -(mu/R^3)*rx - 3*n^2*rx;
g2 = -(mu/R^3-n^2)*ry;
g3 = -rz*(mu*rz/R^3-n^2);
g = [g1;g2;g3];
end
function [Ad, Bd] = CW_dynamics_linear(n, Ts)
    % continuous versions
    Ac = [zeros(3), eye(3); 
        3*n^2, 0, 0, 0, 2*n, 0;
        0, 0, 0, -2*n, 0, 0;
        0, 0, -n^2, 0, 0, 0];
    Bc = [zeros(3);
        eye(3)];
    sys = c2d(ss(Ac, Bc, eye(6), 0), Ts);% discretize
    % discrete vresions
    Ad = sys.A; 
    Bd = sys.B;
end
function [L, f, E, W, Aeq, beq] = prepDockingQP(A, B, N, v_lim, z_lim, p_rel, d_min, v_max, z0)
    nx = 12; % 12 elements in state (6 state vec x 2 agents)
    nu = 6; % 6 elements in control (ux, uy, uz (or u1/2/3 in paper) and 2 agents)
    [H, G] = getHG(A, B, N);
    qBar = 2.25; % tuned form mpcTubeII_QRtuning
    rBar = 1.75; % tuned form mpcTubeII_QRtuning
    Q = eye(nx*N)*qBar; % tune for weight on state
    R = eye(nu*N)*rBar;% tune for weight on control
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
    beq = zeros(0, 1);% active set solver is particular about Aeq, beq
end
function [L, F, E, W, Aeq, beq, H] = prepSKQP(A, B, N, v_lim, z_lim)
    nx = 6; 
    nu = 3;
    [H, G] = getHG(A, B, N);
    Q = eye(nx*N); 
    R = eye(nu*N)*1;
    L = 2*(G'*Q*G + R); 
    F = 2*(G'*Q);
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
    %continuous versions
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

function x_est = kfStep(x0_actual, x0TaskSelection)
x0_ref = x0TaskSelection; % best guess available
num_agents = 8;
n = 0.00113;
t0 = 0;
tf = 5; % consistent with COV implementation, reasonable time-horizon for this mission
tvecPrime = linspace(t0, tf, 1);
M = [eye(3), zeros(3)]; % because position only
W_inv = eye(3) * (2.25^2); % measurement noise
all_norm_pos_errors = cell(num_agents, 1); % will append respective errors for each agent
all_final_pos = cell(num_agents, 1); % will append final positions for each agent
% Define reference trajectory (init guess)
d = 250; % distance parameter for initial spacing

for agentindex = 1:num_agents
    % chaser spacecraft formation
    offsets = x0_actual - x0_actual(:,1); % shift all init pos so agent 1 is at the origin
    % necessary for the dynamics (and translating to C_1 LVLH in KF)
    Re = 6378137;%earth radius m
    h = 500000;%500km altitude
    Rorb = Re + h;
    % Simulate measurements
    yM = cell(agentindex, 1);
    for k = 1:agentindex
        for i = 1:length(tvecPrime)
            Phi_true = cwSTM(tvecPrime(i), n); % Mapping t0 to ti
            x_target_ti = Phi_true * x0_actual(:,k);
            % Measurement is target at ti minus agent offset
            dy = offsets(2,k);
            theta=dy/Rorb;
            RC1Ck=[cos(theta), sin(theta), 0;-sin(theta),cos(theta),0;0, 0, 1];
            initialOffsetTrue = x0_actual(:,k) - x0_actual(:,1);%start update for propagating offset each step
            currentOffsetTrue = Phi_true*initialOffsetTrue;
            dyTrue = currentOffsetTrue(2);
            thetaTrue = dyTrue/Rorb;
            RC1CkTrue = [cos(thetaTrue), sin(thetaTrue), 0; -sin(thetaTrue), cos(thetaTrue), 0; 0, 0, 1];
            yM{k}(:,i) = RC1CkTrue * (M * x_target_ti - currentOffsetTrue(1:3)) + 0.15*randn(3,1);% end update propagaating each step, use DCM to translate to LVLH on C1
        end
    end

    % Init LKF
    x_plus = x0_ref; % x_{i-1}^+
    P_plus = eye(6) * 1e8; % initial covariance  (high init guess)
    xhist = zeros(6, length(tvecPrime));
    for i = 1:length(tvecPrime)
        % STM Mapping - Phi_i is the STM from t_{i-1} to t_i for the time update
        % (for the measurement G, will need Phi from t_0 to t_i)
        dt = 1; % unit step for simplicity
        Phi_step = cwSTM(dt, n);
        Phi_from_t0 = cwSTM(tvecPrime(i), n);
        
        %Reference state at current time
        xi_ref = Phi_from_t0 * x0_ref';
        x_prev_ref = cwSTM(tvecPrime(max(1,i-1)), n) * x0_ref';

        %TIME update (see p352-354 lec notes)
        % x_i^- = x_i^ref + Phi_i(x_{i-1}^+ - x_{i-1}^ref)
        x_minus = xi_ref + Phi_step * (x_plus' - x_prev_ref);
        Q = eye(6) * 1e-9; % model uncertainty
        P_minus = Phi_step * P_plus * Phi_step' + Q;
        % MEASUREMENT Update
        for k = 1:agentindex
            initialOffsetK = x0_actual(:,k) - x0_actual(:,1);
            currentOffsetK = Phi_from_t0*initialOffsetK;
            dy = offsets(2,k);
            theta = dy/Rorb;
            RC1Ck = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
            % Gi = partial g / partial x at reference
            Gi = RC1Ck*M; % DCM must be used for the dynamics
            % z_i: actual measurement
            zi = yM{k}(:,i);
            % g_i(xi_ref): what the measurement would be if state = reference
            % Account for the agent offset in the measurement model
            gi_ref = RC1Ck*(M * xi_ref) -(RC1Ck* currentOffsetK(1:3)); % DCM must be used for the dynamics
            %K_i = P_i^- Gi' (W_inv + Gi P_i^- Gi')^-1
            Si = Gi * P_minus * Gi' + W_inv;
            Ki = P_minus * Gi' / Si;
            %x_i^+ = x_i^- + Ki(zi - gi_ref - Gi(x_minus - xi_ref))
            innov = zi - gi_ref - Gi * (x_minus - xi_ref);
            x_plus = x_minus + Ki * innov;
            % P_i^+ = (I - Ki Gi) * P_i^-
            P_plus = (eye(6) - Ki * Gi) * P_minus;
            % Plus for this timestep is minus for next time step
            x_minus = x_plus;
            P_minus = P_plus;
        end
        x_est = x_minus;
    end
end
end