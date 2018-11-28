function vx = cross(v)
% Gives the matrix form of the cross product

vx = [0 -v(3) v(2)
      v(3) 0 -v(1)
      -v(2) v(1) 0];