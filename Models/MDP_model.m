classdef MDP_model
    % LINMODEL Class of LTI systems with noise on the transitions. 
    % This model class can be used to grid the continuous dynamics to a
    % finite MDP
    
    properties
        type = "MDP abstract";
        states % B matrix for noise
        outputs % noise mean
        inputs
        dim  % the dimension of the state space
        orig % original model
        P
        labels
        beta
        hx
        zstates
        l
        Partition
        outputmap
             
    end
    
    methods
        function obj = MDP_model(P,hx,varargin)
            %MDP_MODEL Construct an instance of this class
            %   Load MDP model as
            %   MDP = MDP_model(P,hx,states,beta, LTI_orig)
            obj.P = P;
            obj.hx = hx;
            if length(varargin)>=1
                obj.states = varargin{1};
            end
            if length(varargin)>=2
                obj.beta = varargin{2}; %!!!! gridsize is twice the size that you would have with a direct gridding. 
            end
            if length(varargin)>=3
                obj.orig = varargin{3};
                obj.outputs = obj.orig.C* obj.states;
            end
        end


        function xn = f_det(obj,x,u)
            xn = obj.orig.f_det(x,u);
        end
        
        
        function xn = sim(obj,x,u)
            % x to index
            % u to index
            
            xn = obj.sim(x,u);
        end
        
        
    end
end
