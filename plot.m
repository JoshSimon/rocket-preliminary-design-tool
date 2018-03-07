    
    % Center of pressure calculation
    
    centerOfPressureNose = L_nose * 0.5;
    centerOfPressureBody = 
    
    
    
    
    % Plot
      
      % reset the graphs
      cla(handles.axes1, 'reset') 
      cla(handles.axes2, 'reset')
      cla(handles.axes3, 'reset')
      cla(handles.axes4, 'reset')
      cla(handles.axes5, 'reset')
    
      % plot the center of pressure and the contour of the rocket
      
      axes(handles.axes1)

        plot(xcp, 0,'Marker','o','MarkerSize',15,'MarkerEdgeColor','red'); hold on;

        % upper body part
        
        plot(linspace(Ln,Lr,10),linspace(R,R,10)); hold on; %Oberseite des Rumpf
        plot(linspace(Ln,xversteck,10),linspace(-R,-R,10)); hold on; %Untereseite des Rumpf

        % lower body part
        
        plot(linspace(LwOVkAp,LwOVkEp,10),linspace(R,Span+R,10));hold on; %obere Leitwerkvorderkante

        plot(linspace(LwOVkEp,LwOHkEp,10),linspace(Span+R,Span+R,10));hold on; %obere Leitwerkspitze

        plot(linspace(LwOHkAp,LwOHkEp,10),linspace(R,Span+R,10));hold on; %obere Leitwerkhinterkante

        plot(linspace(LwOHkAp,LwOHkAp,10),linspace(R,R,10));hold on; %obere Leitwerkunterkante

        plot(linspace(LwOVkAp,LwOVkEp,10),linspace(-RLwV,-RLwV-SpanVer,10));hold on; %Untere Leitwerkvorderkante

        plot(linspace(LwOVkEp,LwOHkEp,10),linspace(-RLwV-SpanVer,-RLwV-SpanVer,10));hold on; %Untere Leitwerkspitze

        plot(linspace(LwOHkAp,LwOHkEp,10),linspace(-RLwV,-RLwV-SpanVer,10));hold on; %Untere Leitwerkhinterkante

        plot(linspace(LwOVkAp,LwOHkAp,10),linspace(-RLwV,-RLwV,10));hold on; %Untere Leitwerkhinterkante

        plot(linspace(Lr,Lr,10),linspace(R,-RLwV,10));hold on; %R�ckseite Rakete

        plot(linspace(Ln,Ln,10),linspace(R,-R,10));hold on; %Trennlinie Nase-Rakete

        plot(NForm);hold on; %Obere Nase
        plot(-NForm);hold on; %Untere Nase

       axis equal;
        
        xlabel('x [mm]'); %Achsenbeschriftung

        ylabel('y [mm]'); %Achsenbeschriftung
        
        xlim ([0 400]);
        ylim ([-100 100]);
