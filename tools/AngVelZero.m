function [AngVel, isterminal, direction] = AngVelZero(tarray,sarray)

AngVel = [sarray(10); sarray(11); sarray(12)]; % The value that needs to be zero
isterminal = [1;1;1]; % halt integration
direction = [0; 0; 0]; % zero can be approached from either direction

end