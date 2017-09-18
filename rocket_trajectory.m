% Program:
% Rocket_Trajectory_Simulation.m
% Multi-stage rocket dynamics and trajectory simulation.
%
% Description:
% Predicts multi-stage rocket dynamics and trajectories based on the given 
% rocket mass, mass flow, engine thrust, launch parameters, and drag coefficient.
% The drag is related to the ambient air pressure which modelled with the ISA
% standard atmosphere.
%
% Variable List:
% Delta = Time step (s)
% t = Time (s)
% Thrust = Thrust (N)
% Mass = Mass (kg)
% Mass_Flow = Flow of mass per second that comes out of the thruster (kg/s)
% Mass_Rocket_With_Motor = Mass with motor (kg)
% Mass_Rocket_Without_Motor = Mass without motor (kg)
% Theta = Angle (deg)
% C = Drag coefficient
% Rho = Air density (kg/m^3)
% A = Rocket projected area (m^2)
% Gravity = Gravity (m/s^2)
% Launch_Rod_Length = Length of launch rod (m)
% n = Counter
% Fn = Normal force (N)
% Drag = Drag force (N)
% Fx = Sum of forces in the horizontal direction (N)
% Fy = Sum of forces in the vertical direction (N)
% Vx = Velocity in the horizontal direction (m/s)
% Vy = Velocity in the vertical direction (m/s)
% Ax = Acceleration in the horizontal direction (m/s^2)
% Ay = Acceleration in the vertical direction (m/s^2)
% x = Horizontal position (m)
% y = Vertical position (m)
% Distance_x = Horizontal distance travelled (m)
% Distance_y = Vertical travelled (m)
% Distance = Total distance travelled (m)
% Memory_Allocation = Maximum number of time steps expected

clear, clc      % Clear command window and workspace

% Simulation parameters
Delta = 0.1;                  % Time step 
Memory_Allocation = 10000;    % Maximum number of time steps expected

% Preallocate memory for arrays
t = zeros(1, Memory_Allocation);
Thrust = zeros(1, Memory_Allocation);
Mass = zeros(1, Memory_Allocation);
Theta = zeros(1, Memory_Allocation);
Fn = zeros(1, Memory_Allocation);
Drag = zeros(1, Memory_Allocation);
Fx = zeros(1, Memory_Allocation);
Fy = zeros(1, Memory_Allocation);
Ax = zeros(1, Memory_Allocation);
Ay = zeros(1, Memory_Allocation);
Vx = zeros(1, Memory_Allocation);
Vy = zeros(1, Memory_Allocation);
x = zeros(1, Memory_Allocation);
y = zeros(1, Memory_Allocation);
PlotRho = zeros(2, Memory_Allocation); 
Distance_x = zeros(1, Memory_Allocation);
Distance_y = zeros(1, Memory_Allocation);
Distance = zeros(1, Memory_Allocation);
Rho=zeros(1, Memory_Allocation);

% Natural constants
C = 0.4;                                % Drag coefficient
Gravity = 9.81;                         % Gravity (m/s^2)

% Rocket parameters
r = 1.5                                 % Rocket fuselage radius (m)
A = 2*pi*r^2;                           % Rocket projected attack area (m^2)
Launch_Rod_Length = 1;                  % Length of launch rod (m)
Mass_Motor_And_Strucutre_One =  2000;                 % Mass of the first stage rocket motor (kg)
Mass_Motor_And_Structure_Two = 500;                   % Mass of the second stage rocket motor (kg)
Mass_Start = 280000;                    % Start mass of the rocket (kg)
Mass_Fuel_One = 3000;                   % Fuel mass of the first stage (kg)
Mass_Fuel_Two = 4000;                   % Fuel mass of the second stage (kg)
Mass_Structure = 0.05 * Mass_Start;     % Mass of the strucutre and the motors (kg)
Mass_Flow_One = 200                     % Propulsion mass flow of the first stage (kg/s)
Mass_Flow_Two = 20                      % Propulsion mass flow of the second stage (kg/s)
Thrust_One = 280000                     % Sum of thrust of the first stage ( (kg*m)/s^2 )
Thrust_Two = 280                        % Sum of thrust of the second stage ( (kg*m)/s^2 )

% Maneuver parameters
Theta(1) = 89;                  % Initial angle (deg)
DeltaTheta = 0.001;              % Flying angle (deg/s)
Height_Start_Gravity_Turn = 80000 % Height when the vertical flight is stopped and the rocket is tilted for the gravity turn maneuver (m)
Eject_One = false;       
Eject_Two = false;

% Inital parameters
Vx(1) = 0;                      % Initial horizontal speed (m/s)
Vy(1) = 0;                      % Initial vertical speed (m/s)
x(1) = 0;                       % Initial horizontal position (m)
y(1) = 0.1;                     % Initial vertical position (m)
Distance_x(1) = 0;              % Initial horizontal distance travelled (m)
Distance_y(1) = 0;              % Initial vertical distance travelled (m)
Distance(1) = 0;                % Initial distance travelled (m)
Rho(1) = 1.2;                   % Initial air density (kg/m^3)
Mass(1) = Mass_Start;           % Initial rocket mass (kg)
Thrust(1) = Thrust_One;         % Inital rocket thrust ( (kg*m)/s^2 )
n = 1;                          % Initial time step            

clc


while y(n) > 0 && n < Memory_Allocation                  % Run until rocket hits the ground
    
    n = n+1;                    % Increment time step
    t(n)= (n-1)*Delta;          % Elapsed time     
    
    % Determine rocket thrust and mass based on launch phase
    if  Eject_One == false && Eject_Two == false 
      stage = 1;
    endif
   
    if Eject_One == true && Eject_Two == false
      stage = 2;
    endif
    
    switch (stage)
     
     case 1 % start phase, stage 1
        Thrust(n) = Thrust_One;
        Mass(n) = Mass_Start -  Mass_Flow_One * t(n);
        if Mass(n) <= ( Mass(n) - Mass_Fuel_One )
          Eject_One = true;
        endif
     
     case 2 % stage 2
          Thrust(n) = Thrust_Two;                         
          Mass(n) = Mass_Start - ( Mass_Motor_And_Structure_One + Mass_Flow_Two * t(n) );
          if Mass(n) <= ( Mass(n) - Mass_Fuel_One )
            Eject_Two = true;
          endif
     
     otherwise % payload separation
        Thrust(n) = 0;
        Mass(n) =  Mass(n-1) - Mass_Motor_And_Structure_Two;
        
   endswitch  
  
    % Normal force calculations  
    if Distance(n-1) <= Launch_Rod_Length       % Launch rod normal force
        Fn(n) = Mass(n)*Gravity*cosd(Theta(1));
    else
        Fn(n) = 0;                              % No longer on launch rod
    endif
    
    % Drag force calculation
    Drag(n)= 0.5*C*Rho(n-1)*A*(Vx(n-1)^2+Vy(n-1)^2); % Calculate drag force
    
    % Sum of forces calculations 
    Fx(n)= Thrust(n)*cosd(Theta(n-1))-Drag(n)*cosd(Theta(n-1)) - Fn(n)*sind(Theta(n-1));                            % Sum x forces
    Fy(n)= Thrust(n)*sind(Theta(n-1))-(Mass(n)*Gravity)- Drag(n)*sind(Theta(n-1))+Fn(n)*cosd(Theta(n-1));    % Sum y forces
        
    % Acceleration calculation
    Ax(n)= Fx(n)/Mass(n);                       % Net accel in x direction 
    Ay(n)= Fy(n)/Mass(n);                       % Net accel in y direction
	
    % Velocity calculations
    Vx(n)= Vx(n-1)+Ax(n)*Delta;                 % Velocity in x direction
    Vy(n)= Vy(n-1)+Ay(n)*Delta;                 % Velocity in y direction
	
    % Position calculations
    x(n)= x(n-1)+Vx(n)*Delta;                   % Position in x direction
    y(n)= y(n-1)+Vy(n)*Delta;                   % Position in y direction
    
    % Distance calculations    
    Distance_x(n) = Distance_x(n-1)+abs(Vx(n)*Delta);      % Distance in x 
    Distance_y(n) = Distance_y(n-1)+abs(Vy(n)*Delta);      % Distance in y 
    Distance(n) = (Distance_x(n)^2+Distance_y(n)^2)^(1/2); % Total distance
     
    
    % Air pressure
    Rho(n) = 2; 
    disp([Distance_x(n), isa(Distance_x(n),"P")]);
   
    
    % Rocket angle calculation
    Theta(n)= atand(Vy(n)/Vx(n));      % Angle defined by velocity vector
    if Distance_x > Height_Start_Gravity_Turn
      Theta(n) = Theta(n) * DeltaTheta;
    endif
    
end