function P_xbs=Power_BS_ratio(w,cl,ct,a,b,k_b,h1)

 p=(w^2/cl^2-k_b^2)^0.5;
 q=(w^2/ct^2-k_b^2)^0.5;
 nElemExtend = 2*h1/100;
 z=linspace(-h1+nElemExtend/2,h1-nElemExtend/2,100); 
 A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比
 u1=(1i*k_b*cos(p*z)+q*A1*cos(q*z));
 u2=(-p*sin(p*z)-1i*k_b*A1*sin(q*z));
 Txx=(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z)));
 Txz=(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z)));
 G=-1i*w*(u1*conj(Txx.')+u2*conj(Txz.'));
 P_xbs=-0.5*nElemExtend*real(G);

