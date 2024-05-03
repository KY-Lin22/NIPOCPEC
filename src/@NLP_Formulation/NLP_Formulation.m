classdef NLP_Formulation < handle
    % formulate a NLP based on the given OCPEC
    %
    % OCPEC has the form:
    %  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
    %  s.t. Dot{x} = f(x, u, lambda)
    %       lambda \in SOL(K, F(x, u, lambda))
    %       K := {lambda | bl <= lambda <= bu}
    %       G(x, u) >= 0,
    %       C(x, u) = 0,
    %
    % NLP has the form:
    %  min  J(z),
    %  s.t. h(z) = 0,
    %       c(z) >= 0,
    %       g(z, s) >= 0
    
    properties
        z % symbolic variable, includes all the variable to be optimized
        s % symbolic variable, relaxation parameter
        J % symbolic function, cost function 
        h % symbolic function, equality constraint    
        c % symbolic function, inequality constraint without relaxation parameter
        g % symbolic function, inequality constraint with relaxation parameter
      
        Dim % struct, problem dimension record
        FuncObj % structure, function object
    end
    
    %% Constructor method
    methods
        function self = NLP_Formulation(OCPEC)
            %NLP_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            
            %% discretize OCPEC into NLP
            nlp = self.createNLP(OCPEC);
            self.z = nlp.z;   
            self.s = nlp.s;
            self.J = nlp.J;    
            self.h = nlp.h;  
            self.c = nlp.c;
            self.g = nlp.g;    
            self.Dim = nlp.Dim;        
            
            %% create function object
            self.FuncObj = self.createFuncObj(nlp);            
 
            %% display NLP information
            disp('*----------------------------------- NLP Information ------------------------------------*')
            disp('1. Problem Size')
            disp(['number of decision variable (z): ....................... ', num2str(self.Dim.z)])
            disp(['number of equality constraint (h): ..................... ', num2str(self.Dim.h)])
            disp(['number of inequality constraint without paramer (c): ... ', num2str(self.Dim.c)])
            disp(['number of inequality constraint with param (g): ........ ', num2str(self.Dim.g)])

        end
    end
    %% Other method
    methods
        nlp = createNLP(self, OCPEC)

        [constraint_without_param, constraint_with_param] = createRelaxedEquilibriumConstraint(self, OCPEC)
      
        FuncObj = createFuncObj(self, nlp)
    end
end

