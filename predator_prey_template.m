function predator_prey_template
   %% NOTE: DO NOT SUBMIT THIS SCRIPT
   %  Submit only the compute_f_groupname() function below, in a
   %  separate file, which must be named compute_f_groupname.m
   %  (try to think of a unique groupname...)
   %%
   close all
   g = 9.81;
   mr = 100; % Mass of predator, in kg
   my = 10.; % Mass of prey, in kg
   Frmax = 1.3*mr*g; % Max force on predator, in Newtons
   Fymax = 1.4*my*g; % Max force on prey, in Newtons
   c = 0.2; % Drag coeft, in N s/m
   initial_w = [150,1000,0,1000,0,0,0,0]; % Initial position/velocity
   force_table_predator = rand(51,2)-0.5;
   force_table_prey = rand(51,2)-0.5;
   options = odeset('Events',@event,'RelTol',0.001);
   [time_vals,sol_vals] = ode45(@(t,w) eom(t,w,mr,my,Frmax,Fymax,c,force_table_predator,force_table_prey), ...
       [0:1:250],initial_w,options);
   animate_projectiles(time_vals,sol_vals);
   pr = sol_vals(:,1:2); py = sol_vals(:,3:4);
   dist = sqrt((pr(:,1) - pr(:,1)).^2 +(pr(:,2) - py(:,2)).^2);
   figure
   plot (time_vals, dist)
   
   prs = [transpose(sol_vals(:,1));transpose(sol_vals(:,2))];
   pys = [transpose(sol_vals(:,3));transpose(sol_vals(:,4))];
   vrs = [transpose(sol_vals(:,5));transpose(sol_vals(:,6))];
   vys = [transpose(sol_vals(:,7));transpose(sol_vals(:,8))];
   
   Frs = zeros(251,1);
   for i=1:251
       Frs(i) = norm(compute_f_groupname(time_vals(i),Frmax,Fymax,1,prs(:,i),vrs(:,i),pys(:,i),vys(:,i)));
   end
   
   hold on
   plot(time_vals, Frs/100);
   
   

   
   
 end
function dwdt = eom(t,w,mr,my,Frmax,Fymax,c,forcetable_r,forcetable_y)
% Extract the position and velocity variables from the vector w
% Note that this assumes the variables are stored in a particular order in w.
        pr=w(1:2); vr=w(5:6); py=w(3:4); vy=w(7:8);
        g = 9.81;
% Compute all the forces on the predator      
        amiapredator = true;
        Fr = compute_f_groupname(t,Frmax,Fymax,amiapredator,pr,vr,py,vy);
% The force table varies between +/- 0.5 so this makes the random force
% vary between +/- 0.2*mr*g     
        Frrand = 0.4*mr*g*compute_random_force(t,forcetable_r);
        Frvisc = -vr*norm(vr)*c;
        Frgrav = -mr*g*[0;1];
        Frtotal = Fr+Frrand+Frvisc+Frgrav;
% Now do the forces on the prey        
        amiapredator = false;
        Fy = compute_f_groupname(t,Frmax,Fymax,amiapredator,pr,vr,py,vy);
% The force table varies between +/- 0.5 so this makes the random force
% vary between +/- 0.2*mr*g     
        Fyrand = 0.4*my*g*compute_random_force(t,forcetable_y);
        Fyvisc = -vy*norm(vy)*c;
        Fygrav = -my*g*[0;1];
        Fytotal = Fy+Fyrand+Fyvisc+Fygrav;
% Write similar code below to call your compute_f_groupname function to 
% compute the force on the prey, determine the random forces on the prey, 
% and determine the viscous forces on the prey      
         %enter the appropriate code here to compute dwdt ;
         dwdt = [vr;vy;Frtotal/mr;Fytotal/my];
    end
 
   function [event,stop,direction] = event(t,w)
   pr = w(1:2); py = w(3:4);
% Event function to stop calculation when predator catches prey
% Write your code here? For the event variable, use the distance between 
% predator and prey.  You could add other events to detect when predator/prey hit 
% the ground as well.  See the MATLAB manual for how to detect and 
% distinguish between multiple events if you want to do this 
          event = 1;
          stop = 1;
          direction = -1;
    end    

 function F = compute_f_groupname(t,Frmax,Fymax,amiapredator,pr,vr,py,vy)
     function [c, ceq]=nonlinineqcon(v)
         c = v;
         ceq = [];
     end
 
 opts = optimoptions('fmincon', 'Display', 'off');
    persistent lastTr lastVy lastAy;
%% This is the function that should be submitted
%  Please fill in the information below:
 % Test time and place: Enter the time and room for your test here 
% Group members: list the names of your group members here
%%
%  Variable definitions
%  t:   Time
%  Frmax: Max force that can act on the predator
%  Fymax: Max force that can act on the prey
%  amiapredator: Logical variable - if amiapredator is true,
%            the function must compute forces acting on a predator.
%            If false, code must compute forces acting on a prey.
%  pr - 2D vector with current position of predator eg pr = [x_r;y_r]
%  vr - 2D vector with current velocity of predator eg vr= [vx_r;vy_r]
%  py - 2D vector with current position of prey py = [x_prey;y_prey]
%  vy - 2D vector with current velocity of prey py = [vx_prey;vy_prey]
%      NB:pr,vr,py,vy are all column vectors
%  F - 2D vector specifying the force to be applied to the object
%      that you wish to control F = [Fx;Fy]
%      The direction of the force is arbitrary, but if the
%      magnitude you specify exceeds the maximum allowable
%      value its magnitude will be reduced to this value 
%      (without changing direction)

if (amiapredator)
    ay = 0;
    if ~isempty(lastTr)
        if lastTr ~= t
            ay = (vy - lastVy)/(t - lastTr);
        else
            ay = lastAy;
        end
    end
    %{
    disp(t);
    disp(lastTr);
    disp(vy);
    disp(lastVy);
    disp(ay);
    %}
    py_expt = @(time) (1/2) * ay * min(time, 2)^2 + vy * time + py;
    ar_required = @(t_int) (2 * (py_expt(t_int) - (vr * t_int) - pr)) / (t_int ^ 2);
    ar_required_mag = @(t_int) norm(ar_required(t_int));
    fr_required = @(t_int) (ar_required(t_int) - [0;-9.81]) * 100;
    fr_required_mag = @(t_int) norm(fr_required(t_int));
    
    t_int_best = 20;
    for t_int = [0.1 : 0.1 : 10]
        if fr_required_mag(t_int) < (1.3 * 100 * 9.8)
            t_int_best = t_int;
            break;
        end
    end    
    F = fr_required(t_int_best);
    if norm(F) > (1.3 * 100 * 9.8) 
        disp("NOPE! Shortening...");
        F = F ./ norm(F);
    end
    
    min_gnd_dist = 1.5;
    ar_max_y = Frmax/100-9.81-0.2*9.81;
    vr_y = vr(2);
    pr_y = pr(2);
    
    traject_min = -0.5*vr_y^2/ar_max_y + pr_y;
    if traject_min < min_gnd_dist
        F = [0;Frmax];    
    end 
    
    if pr_y <= 0
        disp(['PREDATOR GROUND COLLISION @y=', num2str(pr_y)]);
    end
    
    lastTr = t;
    lastVy = vy;
    lastAy = ay;
    %disp(F);
else
    % Code to compute the force to be applied to the prey
    F = [sin(0.2*t);-1];
    F = Fymax*F/norm(F);
end
end


% Add the random force script and animation scripts from Appendix here 


function F  = compute_random_force(t,force_table)
% Computes value of fluctuating random force at time t, where 0<t<250.
% The variable force_table is a 51x2 matrix of pseudo-random
% numbers between -0.5 and 0.5, computed using
%  force_table = rand(51,2)-0.5;
% NB ? THE FORCE TABLE MUST BE DEFINED OUTSIDE THIS FUNCTION
% If you define it in here it fries the ode45 function

F = [interp1([0:5:250],force_table(:,1),t);...
     interp1([0:5:250],force_table(:,2),t)];

 
end

function animate_projectiles(t,sols)
figure
xmax = max(max(sols(:,3)),max(sols(:,1)));
xmin = min(min(sols(:,3)),min(sols(:,1)));
ymax = max(max(sols(:,4)),max(sols(:,2)));
ymin = min(min(sols(:,4)),min(sols(:,2)));
dx = 0.1*(xmax-xmin)+0.5;
dy = 0.1*(ymax-ymin)+0.5;
 
for i = 1:length(t)
    clf
    plot(sols(1:i,3),sols(1:i,4),'LineWidth',2,'LineStyle',...
                  ':','Color',[0 0 1]);
    ylim([ymin-dy ymax+dy]);
    xlim([xmin-dx xmax+dx]);
    hold on
    plot(sols(1:i,1),sols(1:i,2),'LineWidth',2,'LineStyle',':',...
                   'Color',[1 0 0]);
    plot(sols(i,1),sols(i,2),'ro','MarkerSize',11,'MarkerFaceColor','r');
    plot(sols(i,3),sols(i,4),'ro','MarkerSize',5,'MarkerFaceColor','g');
    pause(0.1);
end

end



