classdef ParameterContainer
    % This class essentially acts as a storage struct for less cluttered
    % parameter passing. By not altering any of the values, struct is
    % passed through reference semantics.
    properties (Constant)
        V0 = 0;
        kappa = 0;
        lambda = 0;
        r0 = 0;
        a = 0;
        Rls = 0;
        als = 0;
        isospin = 0;
        N = 0;
        Z = 0;
        l = 0;
        k = 0;
        xmin = 0;
        xmax = 0;
        xmatch = 0;
        scenario = 0;
    end

    methods (Static)
        function obj = createObj(V0, kappa, lambda, r0, a, Rls, als, ...
                                isospin, N, Z, l, k, xmin, ...
                                xmax, xmatch, scenario)
            % Constructs container with given parameters
            obj.V0 = V0; 
            obj.kappa = kappa;
            obj.lambda = lambda; 
            obj.r0 = r0; 
            obj.a = a; 
            obj.Rls = Rls; 
            obj.als = als; 
            obj.isospin = isospin; 
            obj.N = N; 
            obj.Z = Z; 
            obj.l = l; 
            obj.k = k; 
            obj.xmin = xmin; 
            obj.xmax = xmax; 
            obj.xmatch = xmatch;  
            obj.scenario = scenario;
        end
    end
end