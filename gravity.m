function a_grav = gravity(y)

R = 6371000;            % avg earth radius
y = y + R;              % dist from center of earth

G = 6.67408*10^(-11);   % Gravitational constant
M = 5.974*10^(24);      % Earth Mass


a_grav=G*M/(y^2);

end