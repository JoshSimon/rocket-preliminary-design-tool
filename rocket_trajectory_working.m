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
  Delta = 1;                  % Time step 
  Memory_Allocation = 100000;    % Maximum number of time steps expected

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
  Rho=zeros(1, Memory_Allocation);
  Mass_Flow = zeros(1, Memory_Allocation);
  % Natural constants
  C = 0.4;                                % Drag coefficient
  Gravity = 9.81;                         % Gravity (m/s^2)
  RHO_NaBH4 = 1.0740 *10^3;               % Density of the additive (kg/m^3), 100%          
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
  Mass_Start = Mass_Motor_And_Structure_One... % Start mass of the rocket (kg)
  + Mass_Motor_And_Structure_Two...
  + Mass_Fuel_One...
  + Mass_Oxidizer_One...
  + Mass_Fuel_Two...
  + Mass_Oxidizer_Two...
  + Mass_Payload;
  Mass_Fuel_And_Oxidizer_One = Mass_Fuel_One + Mass_Oxidizer_One;
  Mass_Fuel_And_Oxidizer_Two = Mass_Fuel_Two + Mass_Oxidizer_Two;  
  
  % Rocket motor
  Mass_Flow_One = 11*12;                  % Propulsion mass flow of the first stage (kg/s)
  Mass_Flow_Two = 12;                     % Propulsion mass flow of the second stage (kg/s)
  Thrust_One = 40000*12;                  % Sum of thrust of the first stage ( (kg*m)/s^2 )
  Thrust_Two = 500000.5;                   % Sum of thrust of the second stage ( (kg*m)/s^2 )
  
  % Rocket size
  A = 2*pi*r^2;                           % Rocket projected attack area (m^2)
  Volume_Oxidizer_One = Mass_Oxidizer_One / Rho_H202 ;  % Oxidizer volume of the first stage
  Volume_Fuel_Stage_One = Mass_Fuel_One / Rho_RP1 ;
  Volume_Oxidizer_Two = Mass_Oxidizer_Two / Rho_H202 ;
  Volume_Fuel_Stage_Two = Mass_Fuel_Two / Rho_RP1;  
  Tank_Height_Fuel_One = Volume_Fuel_Stage_One / (pi* Tank_Radius^2);
  Tank_Height_Oxidizer_One = Volume_Oxidizer_One / (pi * Tank_Radius^2);
  Tank_Height_Fuel_Two = Volume_Fuel_Stage_Two / (pi* Tank_Radius^2);
  Tank_Height_Oxidizer_Two = Volume_Oxidizer_Two / (pi * Tank_Radius^2);
  disp('Rocket dimensions \n');
  disp('Tank_Height_Fuel_One');disp(Tank_Height_Fuel_One);
  disp('Tank_Height_Oxidizer_One');disp(Tank_Height_Oxidizer_One);
  disp('Tank_Height_Fuel_Two');disp(Tank_Height_Fuel_Two);
  disp('Tank_Height_Oxidizer_Two');disp(Tank_Height_Oxidizer_Two); 
 
  % Maneuver parameters
  Theta(1) = 89;                          % Initial angle (deg)
  DeltaTheta = 5;                         % Flying angle (deg/s)
  Height_Start_Gravity_Turn = 20000;      % Height when the vertical flight is stopped and the rocket is tilted for the gravity turn maneuver (m)
  Eject_One = false;       
  Eject_Two = false;
  Mission_Success = false;

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
  stage = 1;                      % Initial stage (1 for the first one)
  Seperation_One = false;         % Stage one at start not seperated
  Seperation_Two = false;         % Stage two at start nit seperated

  Mass_Flow(1) = Mass_Flow_One * Delta; % Inital mass flow equals the mass flow of the first stage times the delta
  clc % clean the console
  
  while y(n) > 0  && Mission_Success == false      % Run until rocket hits the ground or mission completed
     
      n = n+1;                    % Increment time step
      t(n)= (n-1)*Delta;          % Elapsed time     
      
    % Determine rocket thrust and mass based on launch phase 
    if (stage == 1)
      %disp('Stage 1 active');
      Thrust(n) = Thrust_One;
      Mass_Flow(n) = Mass_Flow_One * Delta;
      %disp([Thrust(n),Mass_Flow(n)]);
    else
      if(stage == 2)
        %disp('Stage 2 active');
        Thrust(n) = Thrust_Two;
        Mass_Flow(n) = Mass_Flow_Two * Delta; 
        %disp([Thrust(n),Mass_Flow(n)]);
      else
        Thrust(n)=0; 
      endif
    endif

    Mass(n) = Mass(n-1) - Mass_Flow(n);     
    disp('Mass');disp(Mass(n) - Mass_Flow(n));
    
    % Seperation logic
    if Mass(n) < (Mass_Start - Mass_Fuel_And_Oxidizer_One) && Seperation_One == false;     % seperation condition first stage
      disp('Seperation first stage');
      Seperation_One = true;
      stage = 2;
      Mass(n) = Mass(n) - 999;
    endif

    if Mass(n) < (Mass_Start - (Mass_Fuel_And_Oxidizer_One+Mass_Motor_And_Structure_One + Mass_Fuel_And_Oxidizer_Two)) && Seperation_Two == false % seperation condition second stage
      disp('Seperation second stage');
      Seperation_Two = true;
      stage = 3;
      Mass(n) = Mass_Payload;
    endif

    
      % Drag force calculation
      % C - Drag coefficient, related to the velocity?
      % Rho - Air density
      % A - Refernce area
      % (Vx(n-1)^2+Vy(n-1)^2) - absolute velocity squared
      Drag(n)= 0.5*C*Rho(n-1)*A*(Vx(n-1)^2+Vy(n-1)^2); % Calculate drag force
      
      % Sum of forces calculations 
      Fx(n)= Thrust(n)*cosd(Theta(n-1))-Drag(n)*cosd(Theta(n-1));                     % Sum x forces
      Fy(n)= Thrust(n)*sind(Theta(n-1))-(Mass(n)*Gravity)- Drag(n)*sind(Theta(n-1));  % Sum y forces
      %if(stage != 1 && stage != 2) 
      %  Fx(n) = 0;
      %  Fy(n) = 0;
      %endif
    
    
      % Acceleration calculation
      Ax(n)= Fx(n)/Mass(n);                       % Net accel in x direction 
      Ay(n)= Fy(n)/Mass(n);                       % Net accel in y direction

      % Velocity calculations
      Vx(n)= Vx(n-1)+Ax(n)*Delta;                 % Velocity in x direction
      Vy(n)= Vy(n-1)+Ay(n)*Delta;                 % Velocity in y direction

      % Position calculations
      x(n)= x(n-1)+Vx(n)*Delta;                   % Position in x direction
      y(n)= y(n-1)+Vy(n)*Delta;                   % Position in y direction   
      
      % Air pressure calculation with ISA model
      Rho(n) = isa(y(n),"R");   
      
      % Rocket angle calculation
      % ^ Vy
      % |   / 
      % |  / angle in between velocity vector and x-direction
      % | /  
      % |/---> Vx
      if y(n) >= Height_Start_Gravity_Turn  % if y(n) > Height_Start_Gravity_Turn  % if a certain heigth is reached, the rocket is actively tilting
        if stage == 1
          Theta(n) = atand(Vy(n)/Vx(n)) - DeltaTheta*(n/100);    % % Angle defined by velocity vector in degrees    
          %disp('Theta is now turning');disp(atand(Vy(n)/Vx(n)) - DeltaTheta);
        else
          Theta(n) = atand(Vy(n)/Vx(n)) - 2*DeltaTheta*(n/100);
        endif
      else
        Theta(n) = 89; 
              
      endif
      disp('HERE');disp(Theta(n));disp(n);
      
      % set the mission goals and end the mission if archieved
      if y(n) > 500000
        Mission_Success = true;
        %disp('Mission success');
      endif

       %disp('Fy');disp(Fy(n));
      
  end


  figure('units','normalized','outerposition',[0 0 1 1]) % Maximize plot window

  % Figure 2
  subplot(3,2,2)
  plot(t(1:n),y(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Time (s)'});
  ylabel({'Height (m)'});
  title({'Trajectory'});
  ylim([0 500500]);
  
  % Figure 6
  subplot(3,2,6)
  plot(t(1:n),Mass(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Time (s)'});
  ylabel({'Mass (kg)'});
  title({'Rocket mass'});
  
  % Figure 3
  subplot(3,2,3)
  plot(y(1:n),Vy(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Height (m)'});
  ylabel({'Vy (m/s)'});
  title({'Horizontal Velocity'});

  % Figure 5
  subplot(3,2,5)
  plot(y(1:n),Vx(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Height (m)'});
  ylabel({'Vy (m/s)'});
  title({'Vertical Velocity'});


  % Figure 1
  subplot(3,2,1)
  plot(y(1:n),Thrust(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Height (m)'});
  ylabel({'Thrust (N)'});
  title({'Thrust' });


  % Figure 4
  subplot(3,2,4)
  plot(t(1:n),Theta(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Time (s)'});
  ylabel({'Theta (Deg)'});
  title({'Theta'});

  figure('units','normalized','outerposition',[0 0 1 1]) % Maximize plot window
  
  
  subplot(2,1,1)
  plot(y(1:n),Drag(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlabel({'Height (m)'});
  ylabel({'Drag (N)'});
  title({'Drag Force'});

  subplot(2,1,2)
  plot(y(1:n),Rho(1:n), "linewidth", 5);
  set(gca, "linewidth", 4, "fontsize", 30); 
  xlim([0 50000]);
  xlabel({'Height (m)'});
  ylabel({'Density (kg/m^3)'});
  title({'Air density'});
