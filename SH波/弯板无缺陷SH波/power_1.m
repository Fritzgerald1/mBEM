function [p_m,p_p,p_m_p]=power_1(arrRes,Ind,Miu,Rho,H,Omega,n,ElemLen_m_power,...
ElemLen_p_power,MidPoint_m_power,MidPoint_p_power,NormalVector_m_power,NormalVector_p_power,r)

cT=(Miu/Rho)^0.5;
NormalVector_m_power = -NormalVector_m_power;
X1=MidPoint_m_power(1,1);
X2=MidPoint_m_power(:,2)-r-1;
X3=MidPoint_p_power(1,1);
X4=MidPoint_p_power(:,2)-r-1;
[nPt,temp]=size(MidPoint_m_power);      
u_Plus = zeros(nPt,n);
u_Minus = zeros(nPt,n);
t_Plus = zeros(nPt,n);
t_Minus = zeros(nPt,n);
u_minus_plus = zeros(nPt,n);
t_minus_plus = zeros(nPt,n);

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
        t_Plus(:,loop4)=u_Stress_Plus_x1.*NormalVector_p_power(:,1)+u_Stress_Plus_x2.*NormalVector_p_power(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*sin(beta*X2)*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=Miu*beta*cos(beta*X2)*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector_m_power(:,1)+u_Stress_Minua_x2.*NormalVector_m_power(:,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*sin(beta*X2)*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=Miu*beta*cos(beta*X2)*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1.*NormalVector_m_power(:,1)+u_Stress_minus_Plus_x2.*NormalVector_m_power(:,2);

    else
        u_Plus(:,loop4)=cos(beta*X4)*exp(1i*kesi*X3);
        u_Minus(:,loop4)=cos(beta*X2)*exp(-1i*kesi*X1);
        u_minus_plus(:,loop4)=cos(beta*X2)*exp(1i*kesi*X1);

        u_Stress_Plus_x1=Miu*1i*kesi*cos(beta*X4)*exp(1i*kesi*X3);
        u_Stress_Plus_x2=-Miu*beta*sin(beta*X4)*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*NormalVector_p_power(:,1)+u_Stress_Plus_x2.*NormalVector_p_power(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*cos(beta*X2)*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=-Miu*beta*sin(beta*X2)*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector_m_power(:,1)+u_Stress_Minua_x2.*NormalVector_m_power(:,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*cos(beta*X2)*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=-Miu*beta*sin(beta*X2)*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1.*NormalVector_m_power(:,1)+u_Stress_minus_Plus_x2.*NormalVector_m_power(:,2);

    end

end

%%能量守恒

R_m=arrRes(Ind(6)+1:Ind(6)+n,1);
R_p=arrRes(Ind(6)+1+n:Ind(6)+2*n,1);

v_m=(-u_Minus*1i*Omega)*R_m;
v_p=(-u_Plus*1i*Omega)*R_p;
v_m_p=-u_minus_plus(:,1)*1i*Omega;

t_m=t_Minus*R_m;
t_p=t_Plus*R_p;
t_m_p=t_minus_plus(:,1);

p_m=(v_m.'*conj(t_m)+conj(v_m.')*t_m)*ElemLen_m_power/4;
p_p=(v_p.'*conj(t_p)+conj(v_p.')*t_p)*ElemLen_p_power/4;
p=p_m+p_p;
p_m_p=(v_m_p.'*conj(t_m_p)+conj(v_m_p.')*t_m_p)*ElemLen_m_power/4;














