function P_xa=Power_IA_ratio(w,cl,ct,a,b,k,h1)

 p=(w^2/cl^2-k^2)^0.5;
 q=(w^2/ct^2-k^2)^0.5;
 nElemExtend = 2*h1/100;
 z=linspace(-h1+nElemExtend/2,h1-nElemExtend/2,100); 
 A2=-2*1i*k*p*cos(p*h1)/((k^2-q^2)*cos(q*h1));%反对称模式振幅比
 u1=(1i*k*sin(p*z)-q*A2*sin(q*z));%%%%%
 u2=(p*cos(p*z)-1i*k*A2*cos(q*z));
 Txx=(a*(-(k^2+p^2)*sin(p*z))+2*b*(-k^2*sin(p*z)-1i*k*q*A2*sin(q*z)));
 Txz=(b*(2*1i*k*p*cos(p*z)+(k^2-q^2)*A2*cos(q*z)));
 G=-1i*w*(u1*conj(Txx.')+u2*conj(Txz.'));
 P_xa=-0.5*nElemExtend*real(G);

