function [Y_k, Info] = LineSearchMerit(self, beta, Y, p, dY)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;

Y_Node = cumsum([self.NLP.Dim.h, self.NLP.Dim.c, self.NLP.Dim.g, self.NLP.Dim.z]);
% load parameter
stepSize_min = self.Option.LineSearch.stepSize_Min;
stepSize_decayRate = self.Option.LineSearch.stepSize_DecayRate;
rho = self.Option.LineSearch.rho;
nu_D = self.Option.LineSearch.nu_D;
timeStep = self.OCPEC.timeStep;

%% some quantities at current iterate Y
% directional derivative of cost 
J_grad_times_dz = full(self.FuncObj.J_grad_times_dz(Y, dY));
% constraint violation M (L1 norm scaled by time step)
% - L1 norm follows IPOPT, and also the cost is the sum of stage cost
% - as a constraint measure, it need to be scaled by time step to consistent with the cost that has been scaled
M = timeStep * norm(full(self.FuncObj.M(Y, p)), 1);
% penalty parameter
beta_Trial = J_grad_times_dz/((1 - rho) * M);
if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end
% merit and its directional derivative
z  = Y(Y_Node(3) + 1 : Y_Node(4), 1);
merit = full(self.FuncObj.J(z)) + beta_k * M;
merit_DD = J_grad_times_dz - beta_k * M;

%% backtracking line search
stepSize_init = 1;
while true
     %% Step 1: estimate trial stepsize, iterate, and merit
     stepSize_trial = max([stepSize_init, stepSize_min]);
     Y_trial = Y + stepSize_trial * dY;     
     z_trial = Y_trial(Y_Node(3) + 1 : Y_Node(4), 1);
     M_trial = timeStep * norm(full(self.FuncObj.M(Y_trial, p)), 1);
     merit_trial = full(self.FuncObj.J(z_trial)) + beta_k * M_trial;
     %% Step 2: check sufficient decrease condition
     if merit_trial <= merit + stepSize_trial * nu_D * merit_DD
         status = 1;
         break
     end
     %% Step 3: checking min stepsize
     if stepSize_trial == stepSize_min
         % fails to satisfy sufficient decrease condition with min stepsize, break backtracking linesearch
         status = 0;
         break
     else
         % estimate a smaller stepsize
         stepSize_init = stepSize_decayRate * stepSize_init;
     end
end

%% organize output
timeElapsed = toc(timeStart);

Info.status = status;
Info.time = timeElapsed;
switch status
    case 0
        % fail, return the previous one
        Y_k = Y;           
    case 1
        % success, return the new iterate
        Y_k = Y_trial;
        Info.beta = beta_k;
        Info.stepSize = stepSize_trial;
        Info.merit = [merit, merit_trial];
end

end