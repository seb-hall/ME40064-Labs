classdef LayerProperties
    properties
        x       double  % min x coordinate for this layer
        D       double  % diffusion coefficient
        beta    double  % extra-vascular diffusivity
        gamma   double  % drug degradation rate
    end

    methods
        function obj = LayerProperties(x, D, beta, gamma)
            obj.x = x;
            obj.D = D;
            obj.beta = beta;
            obj.gamma = gamma;
        end
    end
end