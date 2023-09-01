function P_xba=Power_BA_ratio(f,k_b,h1)
                
 w=2*pi*f;
 cl=5940;               
 ct=3200;
 b=7800*ct^2;
 a=7800*cl^2-2*b;
 p=(w^2/cl^2-k_b^2)^0.5;
 q=(w^2/ct^2-k_b^2)^0.5;
 nElemExtend = 2*h1/100;
 z=linspace(-h1+nElemExtend/2,h1-nElemExtend/2,100); 
 A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比
 u1=(1i*k_b*sin(p*z)-q*A2*sin(q*z));%%%%%原始值没有负号
 u2=(p*cos(p*z)-1i*k_b*A2*cos(q*z));
 Txx=(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-i*k_b*q*A2*sin(q*z)));
 Txz=(b*(2*i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z)));
 G=-1i*w*(u1*conj(Txx.')+u2*conj(Txz.'));
 P_xba=-0.5*nElemExtend*real(G);

