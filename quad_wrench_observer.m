% Reference: https://hal.science/hal-01137966/document
% A nonlinear force observer for quadrotors and
%  application to physical interactive tasks
% Assumptions:
% - RPY given in intrinsic euler
% - Coriolis forces negligible (same assumption made in paper
%   implicitly by neglecting drag forces)

% States
% xb = [x; y; z]; (NED)
% eta = [roll; pitch; yaw];
% eta = [phi; theta; psi];
% q = [xb'; eta']

clear; close all; rng(1);

%% Simulation settings
tf = 50;
dt = 0.001;
t = 0:dt:tf;
n = length(t);

% Measurement noise standard deviations
position_std = 0.15e-3; % m
orientation_std = 0.15e-3; % rad
velocity_std = 1.2e-3; % m/s
rates_std = 1.2e-3; % rad/s
q_noise_matrix = diag([position_std * ones(1,3), orientation_std * ones(1,3)]);
q_dot_noise_matrix = diag([velocity_std * ones(1,3), rates_std * ones(1,3)]);

%% System Dynamics
m = 1;
g0 = 9.81;
Jb = diag([0.13, 0.13, 0.22]);
e3 = [0;0;1];
c = @(x) cos(x);
s = @(x) sin(x);

% Body frame expressed wrt world frame Rwb
Rx = @(x) [1 0 0; 0 c(x) -s(x); 0 s(x) c(x)];
Ry = @(x) [c(x) 0 s(x); 0 1 0; -s(x) 0 c(x)];
Rz = @(x) [c(x) -s(x) 0; s(x) c(x) 0; 0 0 1];
Rwb = @(q) Rx(q(4)) * Ry(q(5)) * Rz(q(6));

% Tranformation matrix
T = @(q) [1 0 -s(q(5)); 
         0 c(q(5)) s(q(4))*c(q(5)); 
         0 -s(q(4)) c(q(4))*c(q(5))];

% Rotational inertia expressed in world frame
Jw = @(q) T(q)'*Jb*T(q);

% Mass matrix
B = @(q) blkdiag(m*eye(3), Jw(q));

% Gravity force
g = [-m*g0*e3; zeros(3,1)];

% Input gain (6x4)
G = @(q) [Rwb(q)*e3 zeros(3,3); 
          zeros(3,1) eye(3)];

%% Observer Dynamics 
% observer gains
c_force_obs_gain = 10;
c_torque_obs_gain = 10.0;
p = @(q_dot, c_obs) c_obs * q_dot;
L = @(q, c_obs) c_obs * eye(6) / inv(B(q));

%% Controller parameters
force_z_desired = 0.0; % Down, N

%% Simple estimator design
tc = 0.2/5; % s
alpha = dt/(dt + tc);
tc_estimate = 0.4/5;

%% Simulate
q = zeros(6, n);
q_dot = zeros(6, n);
q_measured = zeros(6,n);
q_dot_measured = zeros(6,n);
f_x_external = 0.5 + 0.5*square((t-5)*2*pi/30, 5/30*100);
f_z_external = -0.5 - 0.5*square((t-15)*2*pi/20, 5/20*100);
tau_x_external = 0.5 * (0.5 + 0.5*square((t-25)*2*pi/10, 5/10*100));
tau_x_external(1:int64(25/dt)) = 0;
tau_x_external(int64(40/dt):end) = 0;
wrench_external = zeros(6,n);
wrench_external(1,:) = f_x_external;
wrench_external(3,:) = f_z_external;
wrench_external(4,:) = tau_x_external;
u_control = zeros(4,n);
psi_force = zeros(6,n); % [force_external; torque_external]
wrench_force_estimate = zeros(6,n);
psi_torque = zeros(6,n); % [force_external; torque_external]
wrench_torque_estimate = zeros(6,n);
q_ddot_estimate = zeros(6,n);
wrench_simple_estimate = zeros(6,n);
wrench_simple_estimate_filtered = zeros(6,n);
for i = 2:n
    % Compute control effort
    thrust = (force_z_desired - wrench_force_estimate(3,i-1) - m*g0);
    torque_x = -tau_x_external(i-1);
    u_control(:,i-1) = [thrust; torque_x; zeros(2,1)];

    % Integrate states
    q_ddot = B(q(:,i-1)) \ (wrench_external(:,i-1) - g + G(q(:,i-1)) * u_control(:,i-1));
    q_dot(:,i) = q_dot(:,i-1) + dt * (q_ddot);
    q(:,i) = q(:, i-1) + dt * q_dot(:, i);

    % Measurements
    q_measured(:,i) = q(:,i) + q_noise_matrix * randn(6,1);
    q_dot_measured(:,i) = q_dot(:,i) + q_dot_noise_matrix * randn(6,1);

    % Step force observer
    psi_dot_force = L(q_measured(:,i), c_force_obs_gain) * (-psi_force(:,i-1) + g ...
        - G(q_measured(:,i))*u_control(:,i-1) - p(q_dot_measured(:,i), c_force_obs_gain));
    psi_force(:,i) = psi_force(:,i-1) + dt * psi_dot_force;
    wrench_force_estimate(:,i) = psi_force(:,i) + p(q_dot_measured(:,i), c_force_obs_gain);

    % Step torque observer
    psi_dot_torque = L(q_measured(:,i), c_torque_obs_gain) * (-psi_torque(:,i-1) + g ...
        - G(q_measured(:,i))*u_control(:,i-1) - p(q_dot_measured(:,i), c_torque_obs_gain));
    psi_torque(:,i) = psi_torque(:,i-1) + dt * psi_dot_torque;
    wrench_torque_estimate(:,i) = psi_torque(:,i) + p(q_dot_measured(:,i), c_torque_obs_gain);

    % Simple estimator
    q_ddot_finite_diff = (q_dot_measured(:,i) - q_dot_measured(:,i-1)) / dt;
    q_ddot_estimate(:, i) = (1-alpha) * q_ddot_estimate(:, i-1) + alpha * q_ddot_finite_diff;
    wrench_simple_estimate(:, i) = B(q_measured(:,i)) * q_ddot_estimate(:, i) ...
        + g - G(q_measured(:,i)) * u_control(:,i-1);
    wrench_simple_estimate_filtered(:,i) = (1 - dt/(dt+tc_estimate)) ...
        * wrench_simple_estimate_filtered(:,i-1) + dt/(dt+tc_estimate) * wrench_simple_estimate(:,i);
end

%% Plotting
fontsize = 17;
linewidth = 2;
figure(1)

subplot(3,3,1); grid on; hold on;
plot(t, q(1,:), LineWidth=linewidth, color="k")
ylabel("$q_1$ [m]", Interpreter="latex", FontSize=fontsize)
title("Position (x)"); xlim([0,50])

subplot(3,3,2); grid on; hold on;
plot(t, q(2,:), LineWidth=linewidth, color="k")
ylabel("$q_2$ [m]", Interpreter="latex", FontSize=fontsize)
title("Position (y)"); xlim([0,50])

subplot(3,3,3); grid on; hold on;
plot(t, q(3,:), LineWidth=linewidth, color="k")
ylabel("$q_3$ [m]", Interpreter="latex", FontSize=fontsize)
title("Position (z)"); xlim([0,50])

subplot(3,3,4); grid on; hold on;
plot(t, rad2deg(q(4,:)), LineWidth=linewidth, color="k")
ylabel("$\phi$ [deg]", Interpreter="latex", FontSize=fontsize)
title("Orientation (roll)"); xlim([0,50])

subplot(3,3,5); grid on; hold on;
plot(t(1:end-1), -u_control(1,1:end-1), LineWidth=linewidth, color="red")
ylabel("$\rho$ [N]", Interpreter="latex", FontSize=fontsize)
title("Thrust"); xlim([0,50])

subplot(3,3,6); grid on; hold on;
plot(t(1:end-1), u_control(2,1:end-1), LineWidth=linewidth, color="red")
ylabel("$\bar{\tau}_1$ [Nm]", Interpreter="latex", FontSize=fontsize)
title("Torque"); xlim([0,50])

subplot(3,3,7); grid on; hold on;
plot(t, wrench_external(1,:), LineWidth=linewidth, color="k")
plot(t, wrench_force_estimate(1,:), LineWidth=linewidth, color="red")
plot(t, wrench_simple_estimate_filtered(1,:), LineWidth=linewidth, color="blue")
ylabel("$\hat{f}_{e_{x}}$ [N]", Interpreter="latex", FontSize=fontsize)
title("Observer Estimation")
xlabel("Time [s]"); xlim([0,50])

subplot(3,3,8); grid on; hold on;
plot(t, wrench_external(3,:), LineWidth=linewidth, color="k")
plot(t, wrench_force_estimate(3,:), LineWidth=linewidth, color="red")
plot(t, wrench_simple_estimate_filtered(3,:), LineWidth=linewidth, color="blue")
ylabel("$\hat{f}_{e_{z}}$ [N]", Interpreter="latex", FontSize=fontsize)
title("Observer Estimation")
xlabel("Time [s]"); xlim([0,50])

subplot(3,3,9); grid on; hold on;
plot(t, wrench_external(4,:), LineWidth=linewidth, color="k", DisplayName="Exact")
plot(t, wrench_torque_estimate(4,:), LineWidth=linewidth, color="red", DisplayName="Est")
plot(t, wrench_simple_estimate_filtered(4,:), LineWidth=linewidth, color="blue", DisplayName="Simple")
ylabel("$\hat{\tau}_{e_{x}}$ [Nm]", Interpreter="latex", FontSize=fontsize)
title("Observer Estimation")
xlabel("Time [s]"); xlim([0,50])
legend(); 

% Measurement noise
figure(2)
subplot(2,2,1)
plot(t, position_std * randn(1, n), color="blue");
xlim([0,30])
ylabel("[m]")
title("noise in position")

subplot(2,2,2)
plot(t, velocity_std * randn(1, n), color="red");
xlim([0,30])
ylabel("[m/s]")
title("noise in velocity")

subplot(2,2,3)
plot(t, orientation_std * randn(1, n), color="blue");
xlim([0,30])
ylabel("[rd]")
title("noise in orientation")

subplot(2,2,4)
plot(t, rates_std * randn(1, n), color="red");
xlim([0,30])
ylabel("[rd/s]")
title("noise in Euler rates")
