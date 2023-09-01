function P_xbs=Power_BS_ratio(f,k_b,h1)
clc;
format long
%%%%%%%%%%%%%%%实际参数
 syms  z                   
 w=2*pi*f;
 cl=5940;               
 ct=3200;
 b=7800*ct^2;
 a=7800*cl^2-2*b;
 p=(w^2/cl^2-k_b^2)^0.5;
 q=(w^2/ct^2-k_b^2)^0.5;
 A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比
 u1=(1i*k_b*cos(p*z)+q*A1*cos(q*z));
 u2=(-p*sin(p*z)-1i*k_b*A1*sin(q*z));
 Txx=(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z)));
 Txz=(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z)));
 G=u1*conj(Txx)+u2*conj(Txz);
 A=int(G,z,-0.5e-3,0.5e-3);
 P_xbs=-0.5*real(-i*w*A);
 %P_xbs=vpa(p_x);
  
end