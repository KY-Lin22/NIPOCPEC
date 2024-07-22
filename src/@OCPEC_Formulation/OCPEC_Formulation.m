classdef OCPEC_Formulation < handle
    %create an OCPEC with the form of:
    %  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
    %  s.t. Dot{x} = f(x, u, lambda)
    %       lambda \in SOL(K, F(x, u, lambda))
    %       K := {lambda | bl <= lambda <= bu}
    %       G(x, u) >= 0,
    %       C(x, u) = 0,
    
    properties
        timeHorizon % time horizon
        nStages % number of discretized stage
        timeStep % discretization time step
        
        x0 % initial state
        
        x % differentiable state
        u % control input
        lambda % algebraic variable
        
        xMax % x upper bound
        xMin % x lower bound
        uMax % u upper bound
        uMin % u lower bound
        
        L_T % terminal cost
        L_S % stage cost 
        
        f % ODE r.h.s function
        F % VI function
        bl % VI set K: lambda lower bound
        bu % VI set K: lambda upper bound
        VISetType char {mustBeMember(VISetType,{...
            'box_constraint',...
            'nonnegative_orthant'...
            })} = 'box_constraint' % type of VI set K
        G % path inequality constraint (only including bound constraint for x and u, other types of constraints should be transferred into equality constraint)
        C % path equality constraint    
        
        Dim % variable dimemsion record
        FuncObj % CasADi function object        
    end
    
    %% Constructor method
    methods
        function self = OCPEC_Formulation(...
                timeHorizon, nStages, timeStep,...
                x0, ...
                x, u, lambda,...
                xMax, xMin, uMax, uMin,...
                L_T, L_S,...
                f, F, bl, bu, VISetType,...
                G, C)
            %OCPEC_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            disp('creating OCPEC...')
            %% formulate OCPEC            
            % time parameter
            self.timeHorizon = timeHorizon;
            self.nStages = nStages;
            self.timeStep = timeStep;
            % initial state
            self.x0 = x0;
            % variable and their bounds
            self.x = x;
            self.u = u;            
            self.lambda = lambda;
            self.xMax = xMax;
            self.xMin = xMin;
            self.uMax = uMax;
            self.uMin = uMin;            
            % cost function
            self.L_T = L_T;
            self.L_S = L_S;  
            % DVI
            self.f = f;
            self.F = F;     
            self.bl = bl;
            self.bu = bu;
            self.VISetType = VISetType;                       
            % inequality and equality path constraint
            self.G = G;
            self.C = C;
            % dim record
            self.Dim = struct(...
                'x', size(x, 1), 'u', size(u, 1), 'lambda', size(lambda, 1),...
                'G', size(G, 1), 'C', size(C, 1));          
            % function object
            self.FuncObj = self.createFuncObj();

            %% display OCPEC informulation
            disp('*---------------------------------- OCPEC Information -----------------------------------*')
            disp('1. time parameter')
            disp(['time horizon: ............................... ', num2str(self.timeHorizon)])
            disp(['discretization stage: ....................... ', num2str(self.nStages)])
            disp(['time step: .................................. ', num2str(self.timeStep)])            
            disp('2. problem size')
            disp(['VI set type: ................................ ', self.VISetType])
            disp(['number of state variable (x): ............... ', num2str(self.Dim.x)])
            disp(['number of control variable (u): ............. ', num2str(self.Dim.u)])
            disp(['number of algebraic variable (lambda): ...... ', num2str(self.Dim.lambda)])
            disp(['number of path inequality constraint (G): ... ', num2str(self.Dim.G)])
            disp(['number of path equality constraint (C): ..... ', num2str(self.Dim.C)])           

            disp('Done!')
        end
    end
    %% Other method
    methods
        FuncObj = createFuncObj(self)
    end
end

