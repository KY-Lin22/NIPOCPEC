function [state, phase, radius] = holf_oscillator_couple(InitState, InitPhase, nStages, timeStep, option)
%% initialization 
phi_d = InitPhase; % desired phase to specify the gait pattern
num_osci = length(phi_d);

state = zeros(2 * num_osci, nStages);
phase = zeros(1 * num_osci, nStages);
radius = zeros(1 * num_osci, nStages);

% setting oscillator parameter
mu = option.mu; % radius of the limit circle
alpha = 50; % determine the speed for the oscilltor to converge to the limit circle (positive)
b = 50; % ensure mu varies smoothly across different half phase (positive)
beta = option.beta; % fraction of stance phase (i.e.,the phase: phi \in [0, pi])
T = option.T; % period of limit circle

% connection matrix
connectMatrix = cell(num_osci, num_osci);
for i = 1 : num_osci
    for j = 1 : num_osci
        if i == j
            R_ij = zeros(2, 2);           
        else
            theta_dij = phi_d(i) - phi_d(j);
            R_ij = [cos(theta_dij), -sin(theta_dij);...
                    sin(theta_dij), cos(theta_dij)];
        end
        connectMatrix{i, j} = R_ij;
    end
end

%%
for n = 1 : nStages
    % evaluate according to the index of oscillator
    for i = 1 : num_osci
        % previous state
        if n == 1
            rho_Prev_i = InitState(1 + (i - 1) * 2 : i * 2, 1);
        else
            rho_Prev_i = state(1 + (i - 1) * 2 : i * 2, n - 1);
        end
        % parameter in dynamics
        rho_Prev_i_norm = norm(rho_Prev_i, 2);
        gamma_Prev_i = pi / (beta * T * (exp(-b * rho_Prev_i(2)) + 1)) +...
                     pi / ((1 - beta) * T * (exp(b * rho_Prev_i(2)) + 1));
        % dynamics
        f_Prev_i = [alpha * (mu^2 - rho_Prev_i_norm^2) * rho_Prev_i(1) + gamma_Prev_i * rho_Prev_i(2);...
                alpha * (mu^2 - rho_Prev_i_norm^2) * rho_Prev_i(2) - gamma_Prev_i * rho_Prev_i(1)];        
        % couple term
        sum_R_ij_rho_i = [0; 0];
        for j = 1 : num_osci
            R_ij = connectMatrix{i, j};
            sum_R_ij_rho_i = sum_R_ij_rho_i + R_ij * rho_Prev_i;
        end
        % current state
        rho_i = rho_Prev_i + (f_Prev_i + sum_R_ij_rho_i) * timeStep;
        state(1 + (i - 1) * 2 : i * 2, n) = rho_i;
        phase(i, n) = atan2(rho_i(2), rho_i(1));
        radius(i, n) = norm(rho_i, 2);
    end    
end


end

