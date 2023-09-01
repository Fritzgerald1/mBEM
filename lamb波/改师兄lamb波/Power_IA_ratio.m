function P_xa=Power_IA_ratio(f,k,h1)
clc;
format long
%%%%%%%%%%%%%%%实际参数
 syms  z                   
 w=2*pi*f;
 cl=5940;               
 ct=3200;
 b=7800*ct^2;
 a=7800*cl^2-2*b;
 p=(w^2/cl^2-k^2)^0.5;
 q=(w^2/ct^2-k^2)^0.5;
 A2=-2*1i*k*p*cos(p*h1)/((k^2-q^2)*cos(q*h1));%反对称模式振幅比
 u1=(1i*k*sin(p*z)-q*A2*sin(q*z));%%%%%
 u2=(p*cos(p*z)-1i*k*A2*cos(q*z));
 Txx=(a*(-(k^2+p^2)*sin(p*z))+2*b*(-k^2*sin(p*z)-i*k*q*A2*sin(q*z)));
 Txz=(b*(2*i*k*p*cos(p*z)+(k^2-q^2)*A2*cos(q*z)));
 G=u1*conj(Txx)+u2*conj(Txz);
 A=int(G,z,-0.5e-3,0.5e-3);
 P_xa=-0.5*real(-i*w*A);
% P_xba=vpa(p_x);
end