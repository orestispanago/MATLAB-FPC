classdef Air
    %AIR Summary of this class goes here
    properties
        temps;
        densities;
        nys;
        ks;
    end

    methods
        function obj = Air()
            obj.temps =[250 300 350 400 450];
            obj.densities = [1.4235 1.1771 1.0085 0.88213 .8770];
            obj.nys = [11.44 15.89 20.92 26.41 32.39]*1E-6;
            obj.ks = [22.3 26.3 30.0 33.8 37.3]*1E-3;
        end

        function rho = rho(obj,t_air)
            % Calculates air density from temperature (kg/m3)
            rho = interp1(obj.temps, obj.densities, t_air,'spline');
        end

        function ny = ny(obj,t_air)
            % Calculates air kinematic viscosity from temperature (m2/s)
            ny = interp1(obj.temps, obj.nys, t_air,'spline');
        end
        
        function k = k(obj,t_air)
            % Calculates air thermal conductivity from temperature (W/m.k)
            k = interp1(obj.temps, obj.ks, t_air,'spline');
        end
    end
end

