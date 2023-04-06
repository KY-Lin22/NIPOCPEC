classdef Biped < DifferentialVariationalInequalities
    %% Biped: a 2D Biped robot with contact bewteen foot and ground
    %   inspired by thowell's paper and source code https://github.com/thowell/motion_planning/blob/main/models/biped.jl
    
    properties
        mass (3, 1) double {mustBePositive, mustBeFinite} = [0.5 + 0.48 * 2; 0.8112; 0.3037]; % 1: torso; 2: thigh; 3: calf
        inertia (3, 1) double {mustBePositive, mustBeFinite} = [0.0029; 0.00709; 0.00398]; % 1: torso; 2: thigh; 3: calf
        linkLength (3, 1) double {mustBePositive, mustBeFinite} = [0.15 + 0.15; 0.2755; 0.308]; % 1: torso; 2: thigh; 3: calf
        linkCenter (3, 1) double {mustBePositive, mustBeFinite} = [0.0342; 0.2176; 0.1445]; % 1: torso; 2: thigh; 3: calf
        mu = 0.1;% friction coefficient in ground 
        jointFriction = 0.1;
    end
    %% Constructor Method for Biped
    methods
        function plant = Biped(timeStep)
            %Biped
            %   Constructor of Class Biped
            %           plant = Biped(timeStep)
            %
            % Argument: 
            %           timeStep: (1 x 1) double          
            %
            % Output: 
            %           plant: a object of Class Biped
            %
            %% construct an object
            tau_Dim = 5 + 2; % control input: 5; tangential friction force 2
                             % tau1: control input for torso
                             % tau2: control input for thigh 1
                             % tau3: control input for calf 1
                             % tau4: control input for thigh 2
                             % tau5: control input for calf 2
                             % tau6: friction force for leg 1
                             % tau7: friction force for leg 2
            x_Dim = 2 * (2 + 5); % configuration dim : 7
                                 % q1: x pos, 
                                 % q2: z pos
                                 % q3: torso angle (rel. to upward vertical)
                                 % q4: thigh 1 angle (rel. to downward vertical)
                                 % q5: calf 1 (rel. to downward vertical)
                                 % q6: thigh 2 (rel. to downward vertical)
                                 % q7: calf 2 (rel. to downward vertical)
            p_Dim = 6; % 2 contact point, each point need 1 p for normal contact force and 2 p for tangential contact velocity
                       % p1: leg 1 normal contact force
                       % p2, p3: two auxiliary variable for vel in leg 1:  p2 - p3 = velT_1 
                       % p4: leg 2 normal contact force
                       % p5, p6: two auxiliary variable for vel in leg 2:  p5 - p6 = velT_2 
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);            
            
            % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [30; 30; 30; 30; 30;...
                                        1000; 1000];
            plant.DynVarLimit.tau_Min = [-30; -30; -30; -30; -30;...
                                       -1000; -1000]; 
            plant.DynVarLimit.x_Max = [1; 1; 1/2*pi; 1/2*pi; 1/2*pi; 1/2*pi; 1/2*pi;...
                                       5; 5; 5; 5; 5; 5; 5]; 
            plant.DynVarLimit.x_Min = [0; 0; -1/2*pi; -1/2*pi; -1/2*pi; -1/2*pi; -1/2*pi;...
                                      -5; -5; -5; -5; -5; -5; -5]; 
                                  
            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end
            end 
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.computeStateEquationMethod = 2;
            plant.setDynamics()             
        end
        
    end
    %% Other Methods for Biped
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % set initial configuration
        q = setInitConfiguration(plant, x, xita_torso, xita_thigh_1, xita_calf_1, xita_thigh_2)      
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)
        
        plotSimuResult(plant, simuTimeStep, InitState, tau, x, p) 
    end
    
end

