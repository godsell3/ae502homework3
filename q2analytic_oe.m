% 
% Function to compute Keplerian orbital elements from provided initial
% conditions using the EOM derived in Question 2.
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
function[M_, w_, RA_] = q2analytic_oe(conds, rot, tspan, steps)

    %initial conditions
    l0 = conds(1);
    g0 = conds(2);
    h0 = conds(3);
    L0 = conds(4);
    G0 = conds(5);
    H0 = conds(6);
    
    %propagate over time
    const_t = ones(steps,1);
    g = g0*const_t;
    L = L0*const_t;
    G = G0*const_t;
    H = H0*const_t;
    h = rot*(tspan') + h0;
    l = 1./(L.^3).*(tspan') + l0;
    
    %convert to Keplerian elements
    M_ = l;
    w_ = g;
    RA_ = h - g;

end