classdef MeshLayer
    properties
        x               double  % min x coordinate for this layer
        density_ratio   double  % density ratio for this layer
        D               double  % diffusion coefficient
        beta            double  % extra-vascular diffusivity
        gamma           double  % drug degradation rate

        element_count   uint64  % number of elements in this layer
        layer_offset    uint64  % starting element index for this layer
    end

    methods
        function obj = MeshLayer(x, D, beta, gamma, density_ratio)
            obj.x = x;
            obj.D = D;
            obj.beta = beta;
            obj.gamma = gamma;
            obj.density_ratio = density_ratio;

            obj.element_count = 0;
        end
    end
end