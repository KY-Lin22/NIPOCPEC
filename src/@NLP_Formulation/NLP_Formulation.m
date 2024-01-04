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
        
        SparsityPattern % struct, record constraint Jacobian and Lagrangian Hessian sparsity pattern       
        Dim % struct, problem dimension record
        FuncObj % structure, function object
    end
    
    %% Constructor method
    methods
        function self = NLP_Formulation(OCPEC)
            %NLP_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            
            %% discretize OCPEC into NLP
            switch OCPEC.VISetType
                case 'box_constraint'
                    nlp = self.createNLP_box_constraint(OCPEC);
                case 'nonnegative_orthant'
                    nlp = self.createNLP_nonnegative_orthant(OCPEC);
            end           
            self.z = nlp.z;   
            self.s = nlp.s;
            self.J = nlp.J;    
            self.h = nlp.h;  
            self.c = nlp.c;
            self.g = nlp.g;
            self.SparsityPattern = struct(...
                'h_grad', nlp.h_grad,...
                'c_grad', nlp.c_grad,...
                'g_grad', nlp.g_grad,...
                'LAG_hessian', nlp.LAG_hessian);          
            self.Dim = nlp.Dim;        
            
            %% create function object
            self.FuncObj = self.createFuncObj(nlp);
            
            %% display sparsity pattern of constraint Jacobian and Lagrangian Hessian
            figure(100)
            spy(self.SparsityPattern.h_grad)
            title('sparsity pattern: equality constraint Jacobian')  
            
            figure(101)
            spy(self.SparsityPattern.c_grad)
            title('sparsity pattern: inequality (without s) constraint Jacobian')  
            
            figure(102)
            spy(self.SparsityPattern.g_grad)
            title('sparsity pattern: inequality (with s) constraint Jacobian') 
            
            figure(103)
            spy(self.SparsityPattern.LAG_hessian)
            title('sparsity pattern: Lagrangian Hessian')         
            
        end
    end
    %% Other method
    methods
        nlp = createNLP_box_constraint(self, OCPEC)
        
        nlp = createNLP_nonnegative_orthant(self, OCPEC)
        
        FuncObj = createFuncObj(self, nlp)
    end
end

