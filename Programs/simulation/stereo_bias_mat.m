function [ b ] = stereo_bias_mat( x )

x = wrapToPi(x);
b = x;
b(x <= pi/2 & x >= -pi/2) = -x(x <= pi/2 & x >= -pi/2);
b(x > pi/2) = -(pi - x(x > pi/2) );
b(x < -pi/2 ) = x( x < -pi/2 ) + pi;

b=b/pi*2;

