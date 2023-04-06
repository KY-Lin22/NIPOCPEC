function setDynamics(plant)
%setDynamics  
%   Detailed explanation goes here

% load parameter
alpha1 = plant.convergeSpeed1;
alpha2 = plant.convergeSpeed2;
mu_rlc = plant.radiusLimitCycle;
T = plant.periodLimitCycle;
beta = plant.dutyFactor;
b = plant.radiusSmoothVaries;
InitPhase = plant.InitPhase;

k_N = plant.k_N;
k_T = plant.k_T;

num_Osci = length(InitPhase);

% connection matrix 
connectMatrix = cell(num_Osci, num_Osci);
for i = 1 : num_Osci
    for j = 1 : num_Osci
        if i == j
            R_ij = zeros(2, 2); 
        else
            InitPhaseDiff_ij = InitPhase(i) - InitPhase(j);
            R_ij = [cos(InitPhaseDiff_ij), -sin(InitPhaseDiff_ij);...
                    sin(InitPhaseDiff_ij), cos(InitPhaseDiff_ij)];
        end
        connectMatrix{i, j} = R_ij;
    end
end

%% set state equation and equilibrium dynamics
f = sym('f', [plant.Dim.x, 1]);
K = sym('K', [plant.Dim.p, 1]);
for i = 1 : num_Osci
    % variable
%     tau_i = plant.tau(i);
    rho_i = plant.x(1 + (i - 1) * 4 : 2 + (i - 1) * 4);
    
    yN_i = plant.p(1 + (i - 1) * 2);
    yT_i = plant.p(2 + (i - 1) * 2); 
    
    % holf oscillator
    r_i = sqrt((rho_i(1))^2 + (rho_i(2))^2);
    gamma_i = pi / (beta * T * (exp(-b * rho_i(2)) + 1)) + ...
        pi / ((1 - beta) * T * (exp(b * rho_i(2)) + 1));
    f_HolfOsc_i = [alpha1 * (mu_rlc^2 - r_i^2) * rho_i(1) + gamma_i * rho_i(2);...
        alpha2 * (mu_rlc^2 - r_i^2) * rho_i(2) - gamma_i * rho_i(1)];
    sum_R_ij_rho_i = [0; 0]; % couple term
    for j = 1 : num_Osci
        R_ij = connectMatrix{i, j};
        sum_R_ij_rho_i = sum_R_ij_rho_i + R_ij * rho_i;
    end
    
    % foot velocity trajectory 
    g_N_st_i = 0;
%     g_N_sw_i = -k_N * sin(pi * rho_i(1)) * (f_HolfOsc_i(1) + sum_R_ij_rho_i(1)) + tau_i;
    g_N_sw_i = -k_N * sin(pi * rho_i(1)) * (f_HolfOsc_i(1) + sum_R_ij_rho_i(1));
    g_T_st_i = 0;
    g_T_sw_i = -k_T * cos(-pi/2 * rho_i(1)) * (f_HolfOsc_i(1) + sum_R_ij_rho_i(1));   
    f_fvt_i = [g_N_st_i * (1 - yN_i) + g_N_sw_i * yN_i;...
        g_T_st_i * (1 - yT_i) + g_T_sw_i * yT_i];
    
    % state equation
    f(1 + (i - 1) * 4 : 4 + (i - 1) * 4) = [f_HolfOsc_i + sum_R_ij_rho_i; f_fvt_i];
    
    % equilibrium constraint
    K(1 + (i - 1) * 2 : 2 + (i - 1) * 2) = [rho_i(2); rho_i(2)];    
end

plant.f = f;
plant.K = K;
l = zeros(plant.Dim.p, 1);
u = ones(plant.Dim.p, 1);
plant.l = l;
plant.u = u;

end

