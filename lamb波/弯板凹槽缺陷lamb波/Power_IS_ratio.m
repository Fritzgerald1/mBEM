function P_xa=Power_IS_ratio(w,cl,ct,a,b,k,h1)

 p=(w^2/cl^2-k^2)^0.5;
 q=(w^2/ct^2-k^2)^0.5;
 nElemExtend = 2*h1/100;
 z=linspace(-h1+nElemExtend/2,h1-nElemExtend/2,100); 
 A1=2*1i*k*p*sin(p*h1)/((k^2-q^2)*sin(q*h1));%对称模式振幅比
 u1=(1i*k*cos(p*z)+q*A1*cos(q*z));
 u2=(-p*sin(p*z)-1i*k*A1*sin(q*z));
 Txx=(-a*(p^2+k^2)*cos(p*z)+2*b*(-k^2*cos(p*z)+1i*k*A1*q*cos(q*z)));
 Txz=(b*(-2*1i*k*p*sin(p*z)+(k^2-q^2)*A1*sin(q*z)));
 G=-1i*w*(u1*conj(Txx.')+u2*conj(Txz.'));
 P_xa=-0.5*nElemExtend*real(G);
