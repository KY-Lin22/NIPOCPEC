classdef Quadruped < DifferentialVariationalInequalities
    %% Quadruped: a 2D quadruped with contact bewteen foot and ground
    %   inspired by thowell's paper and source code: https://github.com/thowell/motion_planning/blob/main/models/quadruped.jl
    
    properties
        mass (3, 1) double {mustBePositive, mustBeFinite} = [4.713 + 4 * 0.696;...
                                                             1.013;...
                                                             0.166]; % 1: torso; 2: thigh; 3: calf
        inertia (3, 1) double {mustBePositive, mustBeFinite} = [0.056579028 + 4 * 0.696 * 0.183^2.0;...
                                                                0.005139339;...
                                                                0.003014022]; % 1: torso; 2: thigh; 3: calf
        linkLength (3, 1) double {mustBePositive, mustBeFinite} = [0.267;...
                                                                   0.2;...
                                                                   0.2]; % 1: torso; 2: thigh; 3: calf
        linkCenter (3, 1) double {mustBePositive, mustBeFinite} = [0.5 * 0.267 + 0.0127;...
                                                                   0.5 * 0.2 - 0.00323;...
                                                                   0.5 * 0.2 - 0.006435]; % 1: torso; 2: thigh; 3: calf
        mu = 0.5;% friction coefficient in ground 
        jointFriction = 0.1;
    end
    %% Constructor Method for Quadruped
    methods
        function plant = Quadruped(timeStep)
            %Quadruped
            %   Constructor of Class Quadruped
            %           plant = Quadruped(timeStep)
            %
            % Argument: 
            %           timeStep: (1 x 1) double          
            %
            % Output: 
            %           plant: a object of Class Quadruped
            %
            
            %% construct an object
            tau_Dim = 4 + 4 + 4; % control input 8 (thigh 4, calf 4); tangential friction force 4
                                 % tau1: control input for thigh 1
                                 % tau2: control input for calf 1
                                 % tau3: control input for thigh 2
                                 % tau4: control input for calf 2
                                 % tau5: control input for thigh 3
                                 % tau6: control input for calf 3
                                 % tau7: control input for thigh 4
                                 % tau8: control input for calf 4
                                 % tau9: friction force for leg 1
                                 % tau10: friction force for leg 2
                                 % tau11: friction force for leg 3
                                 % tau12: friction force for leg 4
            x_Dim = 2 * (2 + 5 + 4); % configuration dim : 11
                                     % q1: x pos, 
                                     % q2: z pos  
                                     % q3: torso
                                     % q4: thigh 1 ang
                                     % q5: calf 1 ang
                                     % q6: thigh 2 ang
                                     % q7: calf 2 ang
                                     % q8: thigh 3 ang
                                     % q9: calf 3 ang
                                     % q10: thigh 4 ang
                                     % q11: calf 4 ang 
                                     % ang all rel.to downward vertical                                     
            p_Dim = 12; % 4 contact point, each point need 1 p for normal contact force and 2 p for tangential contact velocity
                        % p1: leg 1 normal contact force
                        % p2, p3: two auxiliary variable for vel in leg 1:  p2 - p3 = velT_1 
                        % p4: leg 2 normal contact force
                        % p5, p6: two auxiliary variable for vel in leg 2:  p5 - p6 = velT_2    
                        % p7: leg 3 normal contact force
                        % p8, p9: two auxiliary variable for vel in leg 3:  p8 - p9 = velT_3
                        % p10: leg 4 normal contact force
                        % p11, p12: two auxiliary variable for vel in leg 4:  p11 - p12 = velT_4            
            
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [33.5; 33.5; 33.5; 33.5; 33.5; 33.5; 33.5; 33.5;...
                                        1000; 1000; 1000; 1000];           
            
            plant.DynVarLimit.tau_Min = [-33.5; -33.5; -33.5; -33.5; -33.5; -33.5; -33.5; -33.5;...
                                        -1000; -1000; -1000; -1000];           
            
            plant.DynVarLimit.x_Max = [1; 1; 1*pi; 1*pi; 1*pi; 1*pi; 1*pi; 1*pi; 1*pi; 1*pi; 1*pi;...
                                       50; 50; 50; 50; 50; 50; 50; 50; 50; 50; 50];  
            
            plant.DynVarLimit.x_Min = [0; 0; -1*pi; -1*pi; -1*pi; -1*pi; -1*pi; -1*pi; -1*pi; -1*pi; -1*pi;...
                                       -50; -50; -50; -50; -50; -50; -50; -50; -50; -50; -50];      
                                   
            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
            end
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.computeInvM = true;
            plant.setDynamics()                
            
        end
     
    end
    %% Other Methods for Quadruped
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % set init configuration
        q = setInitConfiguration(plant, x, xita)    
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p)
    end
end

