function omega = Weight(th, varargin)
%WEIGHT Calculates the UCASS Scatting Cross Section Weight Function
%   omega is the value of the weight function at azimuthal angle theta.
%   Refer to Smith et al (2019) for details.

p = inputParser;
p.addParameter('lens_angle', 60*pi/180, @isscalar);
p.addParameter('hole_half_angle', 10.7*pi/180, @isscalar);
p.addParameter('mirror_half_angle', 43.8*pi/180, @isscalar);
p.parse(varargin{:});
Lmh = p.Results.lens_angle;         % Mirror lens angle.
Hm = p.Results.mirror_half_angle;   % Mirror half angle.
Hh = p.Results.hole_half_angle;     % Hole half angle.

th = reshape(th, [], 1);
omega = zeros(size(th, 1), 1);
for i=1:size(th, 1)
    th_i = th(i);
    omega(i) = 1/pi * (phi_calc(Hm, th_i) - phi_calc(Hh, th_i));
end

    function phi = phi_calc(H, theta)
    if (-1 < (cos(H)-cos(Lmh)*cos(theta))/(sin(Lmh)*sin(theta))) && ...
            ((cos(H)-cos(Lmh)*cos(theta))/(sin(Lmh)*sin(theta)) < 1)

        phi = acos((cos(H)-cos(Lmh)*cos(theta))/(sin(Lmh)*sin(theta)));

    elseif (cos(H)-cos(Lmh)*cos(theta))/(sin(Lmh)*sin(theta)) >= 1

        phi = 0;

    elseif (cos(H)-cos(Lmh)*cos(theta))/(sin(Lmh)*sin(theta)) <= -1

        phi = pi;

    end
    end
end