% 
% Function to solve for the eccentric anomaly of an orbit, following Curtis
% Algorithm 3.1.
% INPUTS
%  M_ - mean anomaly (rad)
%  e_ - eccentricity (unitless)
% OUTPUTS
%  E_ - eccentric anomaly (rad)
function[E_] = calc_Eanom(M_, e_)
    
    %set an error tolerance
    err = 1.e-8;

    %find starting value for E
    if (M_)<pi
        E_ = M_ + e_/2;
    else
        E_ = M_ - e_/2;
    end

    %iterate Equation 3.14 from Curtis to find E
    ratio = 1;
    while abs(ratio)>err
        ratio = (E_ - e_*sin(E_) - M_)/(1 - e_*cos(E_));
        E_ = E_ - ratio;
    end
end