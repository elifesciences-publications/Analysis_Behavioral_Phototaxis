function [ b ] = stereo_bias_string( x )
b=x;
for n=1:size(x,2)
    if x(n)<=pi/2 && x(n)>=-pi/2
        b(n)=-x(n);
    end
    if x(n)>pi/2
        b(n)=-(pi-x(n));
    end
    if x(n)<-pi/2
        b(n)=(pi+x(n));
    end
end
b=b/pi*2;
end


