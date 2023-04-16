% 
% ODE used to solve for an orbit given a state vector of r and v vectors
% and gravitational parameter. Incorporates the J2 perturbation.
% INPUTS
%  state - condition vector containing conditions at a specific time
%          [l g h L G H]' (units of DU, TU, radians)
%  omega - perturbation term; here, the rotation of the frame (1/TU)
% OUTPUTS
%  endcond - condition vector after taking a time step
%            [l g h L G H] (units of DU, TU, radians)
function[endcond] = q2ode(~, state, omega)
    L = state(4);
    dl = 1/(L^3);
    dg = 0;
    dh = omega;
    dL = 0;
    dG = 0;
    dH = 0;

    endcond = [dl; dg; dh; dL; dG; dH];
end
