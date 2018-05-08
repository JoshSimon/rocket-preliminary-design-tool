function ret =  rocketParams()

 % Rocket tanks
 
    % Stage 1
    Volume_Oxidizer_One = Mass_Oxidizer_One / Rho_H202 ;  % Oxidizer volume of the first stage, calculated from the required mass and the density of the oxidizer
    Volume_Fuel_One = Mass_Fuel_One / Rho_RP1 ;           % Fuel volume of the first stage, calculated ...
    Tank_Height_Fuel_Two = Volume_Fuel_Stage_Two / (pi* Tank_Radius^2);     % Tank height of the fuel tank of the first stage, calculated by the volume and the set radius
    Tank_Height_Oxidizer_Two = Volume_Oxidizer_Two / (pi * Tank_Radius^2);  % Tank height of the oxidizer tank, ...
    
    
    
    Tank_Height_Fuel_Two = Volume_Fuel_Stage_Two / (pi* Tank_Radius^2);
    Tank_Height_Oxidizer_Two = Volume_Oxidizer_Two / (pi * Tank_Radius^2);
    
    
    
    Volume_Oxidizer_Two = Mass_Oxidizer_Two / Rho_H202 ;
    Volume_Fuel_Stage_Two = Mass_Fuel_Two / Rho_RP1;  
    
    Tank_Height_Fuel_One = Volume_Fuel_Stage_One / (pi* Tank_Radius^2);
    Tank_Height_Oxidizer_One = Volume_Oxidizer_One / (pi * Tank_Radius^2);
  
    disp('Rocket dimensions \n');
    disp('Tank_Height_Fuel_One');disp(Tank_Height_Fuel_One);
    disp('Tank_Height_Oxidizer_One');disp(Tank_Height_Oxidizer_One);
    disp('Tank_Height_Fuel_Two');disp(Tank_Height_Fuel_Two);
    disp('Tank_Height_Oxidizer_Two');disp(Tank_Height_Oxidizer_Two); 

        A = 2*pi*r^2;                           % Rocket projected attack area (m^2)

    
  % Rocket parameters
  r = 0.95;                               % Rocket fuselage radius (m)
  A = 2*pi*r^2;                           % Rocket projected attack area (m^2)
  Launch_Rod_Length = 1;                  % Length of launch rod (m)
  Mass_Motor_And_Structure_One = 999;     % Mass of the first stage rocket motor (kg)
  Mass_Motor_And_Structure_Two = 1151;    % Mass of the second stage rocket motor (kg)
  Mass_Start = 19301;                     % Start mass of the rocket (kg)
  Mass_Fuel_One = 17001;                  % Fuel mass of the first stage (kg)
  Mass_Fuel_Two = 1201;                   % Fuel mass of the second stage (kg)
  Mass_Payload = 300;                     % Mass of Payload (kg)

  
  % Natural constants
    RHO_NaBH4 = 1.0740 *10^3;               % Density of the additive (kg/m^3), 100%          i
    RHO_Water = 1000;                       % Density of water (kg/m^3)
    Rho_H202 = (1.0740 * 10^3)*0.95 + 0.05*1000; % Density of Hydrogenperoxide (kg/m^3), with 5% Water
    Rho_RP1 = 1000*0.95 + 0.05*RHO_NaBH4;   % Density of RP-1 (kg/m^3) with 5% NaBH4 as an additive
    
  % Rocket parameters
    r = .75;                                % Rocket fuselage radius (m)
    Tank_Radius = .7;                                 % Rocket tank radius
    Hull_Thickness = 0.05;                  % Thickness of the Hull, length between hull and tank

  % Rocket masses
    Mass_Motor_And_Structure_One = 500; % Mass of the first stage rocket motor (kg)
    Mass_Motor_And_Structure_Two = 400;    % Mass of the second stage rocket motor (kg)            
    Mass_Fuel_One = 2418.306;               % Fuel mass of the first stage fuel (kg)
    Mass_Oxidizer_One = 14993.49;           % Oxidizer mass of the first stage (kg)
    Mass_Additiv_Fuel_One = 2418.306 * 0.05;   % Mass of the katalytic additive, in this case 5% (kg)
    Mass_Fuel_Two = 459.9722;               % Fuel mass of the second stage fuel (kg)
    Mass_Oxidizer_Two = 2851.828;           % Mass of the oxidizer of the second stage (kg)
    Mass_Payload = 300;                     % Mass of Payload (kg)
    Mass_Start = Mass_Motor_And_Structure_One + Mass_Fuel_One + Mass_Oxidizer_One + Mass_Fuel_Two + Mass_Oxidizer_Two + Mass_Payload;
    Mass_Fuel_And_Oxidizer_One = Mass_Fuel_One + Mass_Oxidizer_One;
    Mass_Fuel_And_Oxidizer_Two = Mass_Fuel_Two + Mass_Oxidizer_Two;  

  % Rocket motor
    Mass_Flow_One = 11*12;                  % Propulsion mass flow of the first stage (kg/s)
    Mass_Flow_Two = 12;                     % Propulsion mass flow of the second stage (kg/s)
    Thrust_One = 30000*12;                  % Sum of thrust of the first stage ( (kg*m)/s^2 )
    Thrust_Two = 38000.5;                   % Sum of thrust of the second stage ( (kg*m)/s^2 )

  
  % Maneuver parameters
  DeltaTheta = 7;                         % Flying angle change per incremental step (deg/s)
  Height_Start_Gravity_Turn = 5000;       % Height when the vertical flight is stopped and the rocket is tilted for the gravity turn maneuver (m)
  Coasting_Phase_After_Stage_One = 5;     % Coasting phase after stage one in (s)


endfunction