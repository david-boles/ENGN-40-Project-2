function F = compute_f_group9(t,Frmax,Fymax,amiapredator,pr,vr,py,vy)
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
function d2ydx2=bicubic_start_acceleration(Xs, Ys, Ms, Xe, Ye, Me)
        h = Ms * (Xe - Xs);
        i = Ye - Ys;
        k = Me * (Xe - Xs);
        b = 3*i - 2*h - k;
        d2ydx2 = (2*b)/((Xe - Xs)^2);
end
function a=compute_acceleration_for_target(pS, vS, pE, vE, time)
        Ax = bicubic_start_acceleration(0, pS(1), vS(1), time, pE(1), vE(1));
        Ay = bicubic_start_acceleration(0, pS(2), vS(2), time, pE(2), vE(2));
        a = [Ax; Ay];
end
function y=sigmoid(x, c)
    if x > (5/c)
        y=1;
    elseif x < (-5/c)
        y=0;
    else
        y=exp(x*c)/(1+exp(x*c));
    end
end
function val=transition_between(pos, val_neg, val_pos)
    s = sigmoid(pos, 100);
    val = s*val_pos + (1 - s)*val_neg;
end
function r=attenuate_when_down(a)
    if a(2) < (-0.7 * 100 * 9.8 * 1.3)
        a(2) = (-0.7 * 100 * 9.8 * 1.3);
    end
    r=a;
end

    persistent lastTr lastVy lastAy;
if (amiapredator)

    ay = 0;
    if ~isempty(lastTr)
        if lastTr ~= t
            ay = (vy - lastVy)/(t - lastTr);
        else
            ay = lastAy;
        end
    else
        lastTr = -1;
    end
    %{
    disp(t);
    disp(lastTr);
    disp(vy);
    disp(lastVy);
    disp(ay);
    %}
    

    if sqrt((pr(1) - py(1))^2 +(pr(2) - py(2))^2) > 100
        %py_expt = @(time) (1/2) * ay * min(time, 0)^2 + vy * min(time,20) + py;
        py_expt = @(time) (1/2) * ay * min(time, 5)^2 + vy * min(time,25) + py;
        vy_expt = @(time)  ay * min(time, 5) + vy;
        ar_required = @(t_int) compute_acceleration_for_target(pr, vr, py_expt(t_int), vy_expt(t_int), t_int);
        ar_required_mag = @(t_int) norm(ar_required(t_int));
        fr_required = @(t_int) (ar_required(t_int) - [0;-9.81]) * 100;
        fr_required_mag = @(t_int) norm(fr_required(t_int));
        t_int_best = 250;
        t_search_max = 250;
    else

        py_expt = @(time) (1/2) * ay * min(time, 5)^2 + vy * min(time,15) + py;

        ar_required = @(t_int) attenuate_when_down(2 * (py_expt(t_int) - (vr * t_int) - pr)) / (t_int ^ 2);
        ar_required_mag = @(t_int) norm(ar_required(t_int));
        fr_required = @(t_int) (ar_required(t_int) - [0;-9.81]) * 100;
        fr_required_mag = @(t_int) norm(fr_required(t_int));
        t_int_best = 10;
        t_search_max = 10;
    end
    
    for t_int = [0.1 : 0.1 : t_search_max]
        if fr_required_mag(t_int) < (1.3 * 100 * 9.81)
            t_int_best = t_int;
            break;
        end
    end    
    F = fr_required(t_int_best);
    if norm(F) > (1.3 * 100 * 9.8) 
        %disp("NOPE! Shortening...");

        F = Frmax*F./ norm(F);

    end
    
    
    min_gnd_dist = 1;
    ar_max_y = Frmax/100-9.81-0.2*9.81;
    vr_y = vr(2);
    pr_y = pr(2);
    traject_min = -0.5*vr_y^2/ar_max_y + pr_y;
    pull_up = transition_between(vr_y, [0;Frmax], F);
    F = transition_between(traject_min - min_gnd_dist, pull_up, F);
    if norm(F) > (1.3 * 100 * 9.8)
        F = Frmax*F./ norm(F);
    end
    
    if pr_y <= 0
        %disp(['PREDATOR GROUND COLLISION @y=', num2str(pr_y)]);
    end 
    
    
    if floor(t) > floor(lastTr)
        %disp(t);
    end
    
    if norm(F) > (1.3 * 100 * 9.8) 
        %disp("NOPE! Shortening...");

        F = Frmax*F./ norm(F);

    end
    
    if(abs(t-lastTr) > 0.01)
        lastTr = t;
        lastVy = vy;
        lastAy = ay;
    end
    
else
       % Code to compute the force to be applied to the prey
    % should try to make prey shake up and down really fast 
    
    if py(2)<70
        F=Fymax*[0; 1];
    else
        rh = pr-py;
        rhmag = norm(rh);
        
        vrel = vy-vr;
        vrmag = norm(vr);
    
        % if distance between predator and prey is less than 40m
        if (rhmag<70)
            Fy = -(vr(1)/(vrmag+1.e-08));
            Fx = (vr(2)/(vrmag+1.e-08));
            if (vrel(1)*Fx+vrel(2)*Fy<0)
                Fy=-Fy;
                Fx=-Fx;
            end

       
        % Gathering the forces together    
        F = Fymax*[Fx;Fy]/norm([Fx;Fy]);
        F = F + 2.*Fymax*[1;0]/py(2);
        F = Fymax*F/norm(F);
        else
           % Sufficiently far apart
           
           % if the direction of the displacement apart is the same as
           % the velocity of the predator
           if ((py-pr)/norm(py-pr)) == vr/norm(vr)
                Fy = -(vr(1)/(vrmag+1.e-08));
                Fx = (vr(2)/(vrmag+1.e-08));
                if (vrel(1)*Fx+vrel(2)*Fy<0)
                    Fy=-Fy;
                    Fx=-Fx;
                end
                F = Fymax*[Fx;Fy]/norm([Fx;Fy]);
                F = F + 2.*Fymax*[1;0]/py(2);
                F = Fymax*F/norm(F);
           else
               
             % Predicting the path of the predator
             dt = 1;
             F = Fymax*(-(py+dt*vy+(pr+dt*vr)))/norm(py+dt*vy-(pr+dt*vr)); 
           end
        end
        
    end

end 


%
%if sqrt((pr(1) - py(1))^2 +(pr(2) - py(2))^2) <= 1
   % disp(['GOTCHA']);
%end




end