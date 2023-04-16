% 
% Function to solve for the position of a satellite over time.
% INPUTS
%  a     - semimajor axis                                   (DU)
%  e     - eccentricity                                     (dimensionless)
%  inc   - inclination                                      (degrees)
%  M     - array of mean anomaly values         [steps x 1] (radians)
%  w     - array of argument of periapse values [steps x 1] (radians)
%  RA    - array of right ascension values      [steps x 1] (radians)
%  steps - number of steps in the simulation
% OUTPUTS
%  rx_ - array of x-coordinates [steps x 1] (DU)
%  ry_ - array of y-coordinates [steps x 1] (DU)
%  rz_ - array of z-coordinates [steps x 1] (DU)
% USES calc_Eanom
function[rx_, ry_, rz_] = calc_r(a, e, inc, M, w, RA, steps)

    %convert to position
    for i=1:steps
        %calculate eccentric anomaly
        E = calc_Eanom(M(i), e);

        %calculate true anomaly
        f = vpa(2*atan(sqrt((1+e)/(1-e))*tan(E/2)));
    
        %find theta
        theta = vpa(w(i) + f);
    
        %calculate r
        r = a*(1 - e^2)/(1 + e*cos(f));
    
        %find r vector
        ri = cos(theta)*cos(RA(i)) - cosd(inc)*sin(RA(i))*sin(theta);
        rj = cos(theta)*sin(RA(i)) + cosd(inc)*cos(RA(i))*sin(theta);
        rk = sind(inc)*sin(theta);
    
        rx_(i) = double(vpa(r*ri)); %DU
        ry_(i) = double(vpa(r*rj)); %DU
        rz_(i) = double(vpa(r*rk)); %DU
    end

end