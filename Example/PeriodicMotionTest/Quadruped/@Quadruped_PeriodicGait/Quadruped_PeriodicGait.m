classdef Quadruped_PeriodicGait < DifferentialVariationalInequalities
    %% Quadruped_PeriodicGait: a 2D quadruped with contact bewteen foot and ground (periodic motion)
    %   - the model of quadruped is inspired by thowell's paper and source code: https://github.com/thowell/motion_planning/blob/main/models/quadruped.jl
    %   - the periodic gait is achieved by the proposed DVI-based framework for the periodic motion
    
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
    
    properties
        convergeSpeed1 = 50 % control the speed for the oscillator to converge to the limit cycle (rho_1)
        convergeSpeed2 = 50 % control the speed for the oscillator to converge to the limit cycle (rho_2)
        radiusLimitCycle = 1 % radius of the limit cycle
        periodLimitCycle = 0.3 % period of the limit cycle
        dutyFactor = 0.4 % duty factor determining the fraction of rho_2 > 0 in a whole period
        radiusSmoothVaries = 50 % ensure radiusLimitCycle varies smoothly across different half planes
        InitPhase = [0; 0; pi; pi]; % determine the gait pattern
        k_N = 0.3 % weight parameter in foot velocity function (normal and tangential)
        k_T = 0.2        
    end
    
    %% Constructor Method for Quadruped_PeriodicGait    
    methods
        function plant = Quadruped_PeriodicGait(timeStep)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %% construct an object           
            tau_Dim = 4 * 2 + 4; 
                        % --8 control input(4 leg, each leg has two input: thigh and cart) 
                        % tau1: control input for thigh 1
                        % tau2: control input for calf 1
                        % tau3: control input for thigh 2
                        % tau4: control input for calf 2
                        % tau5: control input for thigh 3
                        % tau6: control input for calf 3
                        % tau7: control input for thigh 4
                        % tau8: control input for calf 4
                        % --4 tangential friction force
                        % tau9: friction force for leg 1
                        % tau10: friction force for leg 2
                        % tau11: friction force for leg 3
                        % tau12: friction force for leg 4
            x_Dim = 2 * (2 + 5 + 4) + 4 * 4; 
                        % --11 configuration dim
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
                        % --4 periodic DVI for foot trajectory, each periodic DVI has 4 state:
                        % x1 - x2: rho, oscillator state variable
                        % x3 - x4: foot position (normal and tangential)
            p_Dim = 4 * 3 + 4 * 2; 
                        % --4 contact point, each point need 1 p for normal contact force and 2 p for tangential contact velocity
                        % p1: leg 1 normal contact force
                        % p2, p3: two auxiliary variable for vel in leg 1:  p2 - p3 = velT_1 
                        % p4: leg 2 normal contact force
                        % p5, p6: two auxiliary variable for vel in leg 2:  p5 - p6 = velT_2    
                        % p7: leg 3 normal contact force
                        % p8, p9: two auxiliary variable for vel in leg 3:  p8 - p9 = velT_3
                        % p10: leg 4 normal contact force
                        % p11, p12: two auxiliary variable for vel in leg 4:  p11 - p12 = velT_4       
                        % --4 periodic DVI, each periodic DVI needs two auxiliary variables for transforming Filippov DI into DVI
                        % p1 - p2: auxiliary variable (normal and tangential contact velocity)
                        
            plant = plant@DifferentialVariationalInequalities(tau_Dim, x_Dim, p_Dim);   
             
             % initialize dynamics variable limit
            plant.DynVarLimit.tau_Max = [33.5; 33.5; 33.5; 33.5; 33.5; 33.5; 33.5; 33.5;...
                                        1000; 1000; 1000; 1000];           
            
            plant.DynVarLimit.tau_Min = [-33.5; -33.5; -33.5; -33.5; -33.5; -33.5; -33.5; -33.5;...
                                        -1000; -1000; -1000; -1000];           
            
            plant.DynVarLimit.x_Max = [1; 1.2; 1*pi; 1*pi*ones(8, 1);...
                                       50 * ones(11, 1);...
                                       100 * ones(4 * 4, 1)];  
            
            plant.DynVarLimit.x_Min = [0; 0.2; -1*pi; -1*pi*ones(8, 1);...
                                       -50 * ones(11, 1);...
                                       -100 * ones(4 * 4, 1)]; 

            % initialize other properties
            if nargin ~= 0
                if ~isempty(timeStep)
                    plant.timeStep = timeStep;
                end          
            end
            %% set symbolic formulation about dynamics
            plant.codeOptimize = true;
            plant.computeStateEquationMethod = 3;
            plant.qDim = 11;
            plant.setDynamics()                              
            
        end
        
    end
    
    %% Other Methods for Quadruped_PeriodicGait
    methods
        % symbolic formulation about dynamics for codeGen
        setDynamics(plant)
        
        % set init configuration
        q = setInitConfiguration(plant, x, xita)     
        
        % Visualization Method
        plotConfiguration(plant, InitState, RefState)
        
        animateTrajectory(plant, simuTimeStep, InitState, tau, x, p)      
        
    end
    
end

