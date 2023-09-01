function [p_m_G,p_p_G,p_m_p_G]=power_Gauss4cs(arrRes,Ind,Miu,Rho,H,Omega,n,ElemLen_m_power,...
ElemLen_p_power,MidPoint_m_power,MidPoint_p_power,NormalVector_m_power,NormalVector_p_power,arrW,arrX)

cT=(Miu/Rho)^0.5;
NormalVector_m_power = -NormalVector_m_power;
X1=MidPoint_m_power(1,1);
y1=MidPoint_m_power(:,2)-9;
X3=MidPoint_p_power(1,1);
y2=MidPoint_p_power(:,2)-9;
R_m=arrRes(Ind(6)+1:Ind(6)+n,1);
R_p=arrRes(Ind(6)+1+n:Ind(6)+2*n,1);
[nPt,temp]=size(MidPoint_m_power);
u_Plus = zeros(4,n);
u_Minus = zeros(4,n);
t_Plus = zeros(4,n);
t_Minus = zeros(4,n);
u_minus_plus = zeros(4,n);
t_minus_plus = zeros(4,n);
p_mG=zeros(1,nPt);
p_pG=zeros(1,nPt);
p_mpG=zeros(1,nPt);
E=eye(nPt,1);
e=[1,1,1,1].';

   
y_m=arrX.*(ElemLen_m_power/2)+y1(1);
y_p=arrX.*(ElemLen_p_power/2)+y2(1);
X2=y_m;
X4=y_p;
for loop4=1:n
    nMode=loop4-1;
    xx=mod(nMode,2);
    beta=nMode*pi/(2*H);
    kesi=(Omega^2/cT^2-beta^2)^0.5;
    if xx==1
        u_Plus(:,loop4)=sin(beta*X4)*exp(1i*kesi*X3);
        u_Minus(:,loop4)=sin(beta*X2)*exp(-1i*kesi*X1);
        u_minus_plus(:,loop4)=sin(beta*X2)*exp(1i*kesi*X1);

        u_Stress_Plus_x1=Miu*1i*kesi*sin(beta*X4)*exp(1i*kesi*X3);
        u_Stress_Plus_x2=Miu*beta*cos(beta*X4)*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1*NormalVector_p_power(1,1)+u_Stress_Plus_x2*NormalVector_p_power(1,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*sin(beta*X2)*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=Miu*beta*cos(beta*X2)*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1*NormalVector_m_power(1,1)+u_Stress_Minua_x2*NormalVector_m_power(1,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*sin(beta*X2)*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=Miu*beta*cos(beta*X2)*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1*NormalVector_m_power(1,1)+u_Stress_minus_Plus_x2*NormalVector_m_power(1,2);

     else
        u_Plus(:,loop4)=cos(beta*X4)*exp(1i*kesi*X3);
        u_Minus(:,loop4)=cos(beta*X2)*exp(-1i*kesi*X1);
        u_minus_plus(:,loop4)=cos(beta*X2)*exp(1i*kesi*X1);

        u_Stress_Plus_x1=Miu*1i*kesi*cos(beta*X4)*exp(1i*kesi*X3);
        u_Stress_Plus_x2=-Miu*beta*sin(beta*X4)*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1*NormalVector_p_power(1,1)+u_Stress_Plus_x2*NormalVector_p_power(1,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*cos(beta*X2)*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=-Miu*beta*sin(beta*X2)*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1*NormalVector_m_power(1,1)+u_Stress_Minua_x2*NormalVector_m_power(1,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*cos(beta*X2)*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=-Miu*beta*sin(beta*X2)*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1*NormalVector_m_power(1,1)+u_Stress_minus_Plus_x2*NormalVector_m_power(1,2);

     end
end
v_m=(-u_Minus*1i*Omega)*R_m;
v_p=(-u_Plus*1i*Omega)*R_p;
v_m_p=-u_minus_plus(:,1)*1i*Omega;

t_m=t_Minus*R_m;
t_p=t_Plus*R_p;
t_m_p=t_minus_plus(:,1);

p_m_g=zeros(1,4);
p_p_g=zeros(1,4);
p_m_p_g=zeros(1,4);
for i=1:4
    p_m=((v_m(i)*conj(t_m(i))+conj(v_m(i))*t_m(i))*arrW(i)/4)*ElemLen_m_power/2;
    p_m_g(i)=p_m;
    p_p=((v_p(i)*conj(t_p(i))+conj(v_p(i))*t_p(i))*arrW(i)/4)*ElemLen_p_power/2;
    p_p_g(i)=p_p;
    p_m_p=((v_m_p(i)*conj(t_m_p(i))+conj(v_m_p(i))*t_m_p(i))*arrW(i)/4)*ElemLen_p_power/2;
    p_m_p_g(i)=p_m_p;
end
p_m_G=p_m_g*e;
p_p_G=p_p_g*e;    
p_m_p_G=p_m_p_g*e;

a-1

