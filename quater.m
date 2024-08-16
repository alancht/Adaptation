function [vector1]=quater(vector,axis,angle)
% Right handed rotation
q=[cos(angle/2); sin(angle/2)*axis(1); sin(angle/2)*axis(2); sin(angle/2)*axis(3)];
Q=[1-2*q(3)^2-2*q(4)^2, 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));...
   2*(q(2)*q(3)+q(1)*q(4)), 1-2*q(2)^2-2*q(4)^2, 2*(q(3)*q(4)-q(1)*q(2));...
   2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*q(2)^2-2*q(3)^2];

vector1 = Q*vector;