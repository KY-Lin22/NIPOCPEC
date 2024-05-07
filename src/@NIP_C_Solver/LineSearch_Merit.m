function [gamma_h_k, gamma_c_k, gamma_g_k, z_k, Info] = LineSearch_Merit(self, beta,...
    gamma_h, gamma_c, gamma_g, z, dgamma_h, dgamma_c, dgamma_g, dz,...
    s, sigma, J, h, PSI_c, PSI_g, J_grad)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;

% load parameter
stepSize_min = self.Option.LineSearch.stepSize_Min;
stepSize_decayRate = self.Option.LineSearch.stepSize_DecayRate;
rho = self.Option.LineSearch.rho;
nu_D = self.Option.LineSearch.nu_D;

%% some quantities at current iterate z
% directional derivative of cost
J_DD = J_grad * dz;
% constraint violation M (L1 norm scaled by time step)
% - L1 norm follows IPOPT, and also the cost is the sum of stage cost
% - as a constraint measure, it need to be scaled by time step to consistent with the cost that has been scaled
M = self.OCPEC.timeStep * norm([h; PSI_c; PSI_g], 1);
% penalty parameter
beta_Trial = J_DD/((1 - rho) * M);
if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end
% merit and its directional derivative 
merit = J + beta_k * M;
merit_DD = J_DD - beta_k * M;

%% backtracking line search
has_found_new_iterate = false;
stepSize_init = 1;

while ~has_found_new_iterate
     %% Step 1: estimate trial stepsize, iterate, cost, infeasibility and merit
     % step size, multiplier and z
     stepSize_trial = max([stepSize_init, stepSize_min]);
     gamma_h_trial = gamma_h + stepSize_trial * dgamma_h;
     gamma_c_trial = gamma_c + stepSize_trial * dgamma_c;
     gamma_g_trial = gamma_g + stepSize_trial * dgamma_g;
     z_trial       = z       + stepSize_trial * dz;
     % cost
     J_trial = full(self.NLP.FuncObj.J(z_trial));
     % constraint
     h_trial = full(self.NLP.FuncObj.h(z_trial));
     c_trial = full(self.NLP.FuncObj.c(z_trial));
     g_trial = full(self.NLP.FuncObj.g(z_trial, s));
     % FB function
     PSI_c_trial = full(self.FuncObj.PSI_c(gamma_c_trial, c_trial, sigma));
     PSI_g_trial = full(self.FuncObj.PSI_g(gamma_g_trial, g_trial, sigma));
     % constraint infeasibility
     M_trial = self.OCPEC.timeStep * norm([h_trial; PSI_c_trial; PSI_g_trial], 1);
     % merit
     merit_trial = J_trial + beta_k * M_trial;
     
     %% Step 2: check sufficient decrease condition
     if merit_trial <= merit + stepSize_trial * nu_D * merit_DD
         has_found_new_iterate = true;
         status = 1;
     end
     
     %% Step 3: checking min stepsize
    if ~has_found_new_iterate
        if stepSize_trial == stepSize_min
            % linesearch fails on the min stepsize, break backtracking linesearch procedure
            status = 0;
            break
        else
            % estimate a smaller stepsize
            stepSize_init = stepSize_decayRate * stepSize_init;
        end
    end      
end

%% organize output
timeElapsed = toc(timeStart);

Info.status = status;
Info.time = timeElapsed;
switch status
    case 0
        % fail, return the previous one
        gamma_h_k = gamma_h;
        gamma_c_k = gamma_c;
        gamma_g_k = gamma_g;
        z_k = z;      
    case 1
        % success, return the new iterate
        gamma_h_k = gamma_h_trial;
        gamma_c_k = gamma_c_trial;
        gamma_g_k = gamma_g_trial;
        z_k = z_trial;
        Info.J = J_trial;
        Info.h = h_trial;
        Info.c = c_trial;
        Info.g = g_trial;
        Info.beta = beta_k;
        Info.stepSize = stepSize_trial;
        Info.merit = [merit, merit_trial];
end

end

