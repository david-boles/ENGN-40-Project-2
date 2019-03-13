function F = compute_f_groupname(t,Frmax,Fymax,amiapredator,pr,vr,py,vy)
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
    % Code to compute the force to be applied to the predator
    dt = 8;
    F = Frmax*(py+dt*vy-(pr+dt*vr))/norm(py+dt*vy-(pr+dt*vr));
else
    % Code to compute the force to be applied to the prey
    F = [sin(0.2*t);1];
    F = Fymax*F/norm(F);
end
end
