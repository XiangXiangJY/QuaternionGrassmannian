function q = QuaternionReal(q_i)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[a,b,c,d]=parts(q_i);
t=q_i;
a=t.w;
b=t.x;
c=t.y;
d=t.z;
    q=A01232A(a,b,c,d);
end