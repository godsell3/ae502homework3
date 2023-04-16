% 
% Function to compute Keplerian orbital elements from provided initial
% conditions using the EOM derived in Question 1.
% INPUTS
%  conds - array of initial conditions [l g h L G H]' (units of DU, TU,
%          radians)
%  rot   - perturbation term (rotation of frame)      (1/TU)
%  tspan - array of time steps [1 x steps]            (TU)
%  steps - number of steps in the simulation
% OUTPUTS
%  M_  - array of mean anomaly values         [steps x 1] (radians)
%  w_  - array of argument of periapse values [steps x 1] (radians)
%  RA_ - array of right ascension values      [steps x 1] (radians)
function[M_, w_, RA_] = q1analytic_oe(conds, rot, tspan, steps)

    %initial conditions
    lp0 = conds(1);
    gp0 = conds(2);
    hp0 = conds(3);
    Lp0 = conds(4);
    Gp0 = conds(5);
    Hp0 = conds(6);
    
    %propagate over time
    const_t = ones(steps,1);
    gp = gp0*const_t;
    hp = hp0*const_t;
    Lp = Lp0*const_t;
    Gp = Gp0*const_t;
    Hp = Hp0*const_t;
    lp = 1./(Lp.^3).*(tspan') + lp0;
    
    %convert to original variables
    G3 = Gp;
    H3 = Hp;
    g3 = gp;
    L3 = 1*const_t;
    l3 = lp./(1+3*rot*L3.^2.*H3);
    h3 = hp - rot*L3.^3.*l3;
    
    %convert to Keplerian elements
    M_ = l3;
    w_ = g3;
    RA_ = h3 - g3;

end
