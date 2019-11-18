classdef params < matlab.mixin.Copyable
    properties (Constant)
        
        % =====================================================================
        %
        %               FIXED VASCULATURE & BOLD PARAMETERS 
        % =====================================================================
        %
        % For a full explanation for each parameter see Aquino, K et al (2011)
        
        % tissue parameters  
        g       = 0.41;                     % gamma s^(-2)
        rho_f   = 1062;                     % kg m^(-3)
        alpha   = 0.31;                     % unitless
        beta    = 1/params.alpha;           % unitless
        c1      = 6e-8;                     % unitless
        c2      = 1e4*(32)^(-params.beta);  % Pa kg^(-beta) m^(3beta)
        tau     = 1;                        % s
        E_0     = 0.4;                      % unitless
        eta     = params.E_0/params.tau;    % s^(-1)
        psi     = 0.0018;                   % mol kg^(-1)
        V_0     = 0.03;                     % unitless
        Eps_0   = params.V_0*params.rho_f;  % kg m^(-3)
        Q_0     = (params.psi*params.Eps_0)/(1 + 1/(params.eta*params.tau));  % mol m^(-3)                      % 

        % magnetic field parameters at 3T, and at TE=30ms
        k1      = 4.2;                      % unitless
        k2      = 1.7;                      % unitless
        k3      = 0.41;                     % unitless

        % Additonal Terms
        cp      = 1e-7;                     % s^(-1) Pa^(-1)
        L       = 3e-3;                     % m
        k_0     = acos(0.8)/(params.L);     % m^(-1)
        Cz      = 1e-3/(params.k_0^(-1)*sin(params.k_0*params.L)); % may change to L/
        
    end
    
    properties
        
        %%% Default parameters but can be changed
        % tissue
%         kappa   = 0.65;                     % s^(-1)
%         w_f     = 0.56;                     % s^(-1)
        
        % simulation
        Nw      = 2^12;                     % Number of w points
        freqMax = 10;                       % Maximum temporal frequency f
        Nk      = 2^10;                     % Number of k points
        spatialFreqMax = 10000;              % Maximum spatial frequency fk
        
    end
        
end
