function [p_m,p_p,p_m_p]=power_1(a,b,cl,ct,w,h1,rr,num_s,num,t,P_input,k_all,R,ElemLen_m_power,...
ElemLen_p_power,MidPoint_m_power,MidPoint_p_power)

%%%%%%%%%对应不同模态时，得到单位幅值Lamb波的位移和traction 
[nPt,temp]=size(MidPoint_m_power);

u_Minus=zeros(2*nPt,num);
u_Plus=zeros(2*nPt,num);
u_Minus_Plus=zeros(2*nPt,num);
t_Minus=zeros(2*nPt,num);
t_Plus=zeros(2*nPt,num);
t_Minus_Plus=zeros(2*nPt,num);
x=MidPoint_m_power(:,1);
z=MidPoint_m_power(:,2)-rr;
%%%%%在左半部分，各模态单位幅值Lamb波传播方向向左
for tt=1:num
    k_b=-k_all(tt);%%取负值
    p=(w^2/cl^2-k_b^2)^0.5;
    q=(w^2/ct^2-k_b^2)^0.5;
    if tt<=num_s
        P_xbs=Power_BS_ratio(w,cl,ct,a,b,-k_b,h1);
        A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比        
        u11(:,tt)=double((P_input/P_xbs)^0.5*(1i*k_b*cos(p*z)+q*A1*cos(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xbs)^0.5*(-p*sin(p*z)-1i*k_b*A1*sin(q*z)).*exp(1i*k_b*x));     
        Txx_tra(:,tt)=double((P_input/P_xbs)^0.5*(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z))).*exp(1i*k_b*x));
        Txz_tra(:,tt)=double((P_input/P_xbs)^0.5*(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z))).*exp(1i*k_b*x));
    elseif tt>num_s%%%%%反对称模态
        P_xba=Power_BA_ratio(w,cl,ct,a,b,-k_b,h1);
        A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比 
        u11(:,tt)=double((P_input/P_xba)^0.5*(1i*k_b*sin(p*z)-q*A2*sin(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xba)^0.5*(p*cos(p*z)-1i*k_b*A2*cos(q*z)).*exp(1i*k_b*x));
        Txx_tra(:,tt)=double((P_input/P_xba)^0.5*(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
        Txz_tra(:,tt)=double((P_input/P_xba)^0.5*(b*(2*1i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z))).*exp(1i*k_b*x));
        Tzz_tra(:,tt)=double((P_input/P_xba)^0.5*(-a*(k_b^2+p^2)*sin(p*z)-2*b*(p^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));

    end
end
ii=1:nPt;
u_Minus(2*ii-1,:)=u11(ii,:);
u_Minus(2*ii,:)=u22(ii,:);    
t_Minus(2*ii-1,:)=Txx_tra(ii,:);
t_Minus(2*ii,:)=Txz_tra(ii,:);
%%%%%在左半部分，各模态单位幅值Lamb波传播方向向右
for tt=1:num
    k_b=k_all(tt);%%取正值
    p=(w^2/cl^2-k_b^2)^0.5;
    q=(w^2/ct^2-k_b^2)^0.5;
    if tt<=num_s
        A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比
        u11(:,tt)=(1i*k_b*cos(p*z)+q*A1*cos(q*z)).*exp(1i*k_b*x);
        u22(:,tt)=(-p*sin(p*z)-1i*k_b*A1*sin(q*z)).*exp(1i*k_b*x);
        Txx_tra(:,tt)=(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z))).*exp(1i*k_b*x);
        Txz_tra(:,tt)=(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z))).*exp(1i*k_b*x);
        Tzz_tra(:,tt)=(-a*(k_b^2+p^2)*cos(p*z)-2*b*(p^2*cos(p*z)+1i*k_b*q*A1*cos(q*z))).*exp(1i*k_b*x);
    elseif tt>num_s%%%%%反对称模态
        A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比
        u11(:,tt)=(1i*k_b*sin(p*z)-q*A2*sin(q*z)).*exp(1i*k_b*x);
        u22(:,tt)=(p*cos(p*z)-1i*k_b*A2*cos(q*z)).*exp(1i*k_b*x);
        Txx_tra(:,tt)=(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x);
        Txz_tra(:,tt)=(b*(2*1i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z))).*exp(1i*k_b*x);
        Tzz_tra(:,tt)=(-a*(k_b^2+p^2)*sin(p*z)-2*b*(p^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x);
    end
end
ii=1:nPt;
u_Minus_Plus(2*ii-1,:)=u11(ii,:);
u_Minus_Plus(2*ii,:)=u22(ii,:);    
t_Minus_Plus(2*ii-1,:)=Txx_tra(ii,:);
t_Minus_Plus(2*ii,:)=Txz_tra(ii,:);

x=MidPoint_p_power(:,1);
z=MidPoint_p_power(:,2)+rr;
%%%%%在右半部分，各模态单位幅值Lamb波传播方向向右
for tt=1:num
    k_b=k_all(tt);%%取正值
    p=(w^2/cl^2-k_b^2)^0.5;
    q=(w^2/ct^2-k_b^2)^0.5;
    if tt<=num_s
        P_xbs=Power_BS_ratio(w,cl,ct,a,b,k_b,h1);
        A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比
        u11(:,tt)=double((P_input/P_xbs)^0.5*(1i*k_b*cos(p*z)+q*A1*cos(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xbs)^0.5*(-p*sin(p*z)-1i*k_b*A1*sin(q*z)).*exp(1i*k_b*x));
        Txx_tra(:,tt)=double((P_input/P_xbs)^0.5*(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z))).*exp(1i*k_b*x));
        Txz_tra(:,tt)=double((P_input/P_xbs)^0.5*(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z))).*exp(1i*k_b*x));
        Tzz_tra(:,tt)=double((P_input/P_xbs)^0.5*(-a*(k_b^2+p^2)*cos(p*z)-2*b*(p^2*cos(p*z)+1i*k_b*q*A1*cos(q*z))).*exp(1i*k_b*x));
    elseif tt>num_s%%%%%反对称模态
        P_xba=Power_BA_ratio(w,cl,ct,a,b,k_b,h1);
        A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比
        u11(:,tt)=double((P_input/P_xba)^0.5*(1i*k_b*sin(p*z)-q*A2*sin(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xba)^0.5*(p*cos(p*z)-1i*k_b*A2*cos(q*z)).*exp(1i*k_b*x));
        Txx_tra(:,tt)=double((P_input/P_xba)^0.5*(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
        Txz_tra(:,tt)=double((P_input/P_xba)^0.5*(b*(2*1i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z))).*exp(1i*k_b*x));
        Tzz_tra(:,tt)=double((P_input/P_xba)^0.5*(-a*(k_b^2+p^2)*sin(p*z)-2*b*(p^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
    end
end
u_Plus(2*ii-1,:)=u11(ii,:);
u_Plus(2*ii,:)=u22(ii,:);    
t_Plus(2*ii-1,:)=Txx_tra(ii,:);
t_Plus(2*ii,:)=Txz_tra(ii,:);
%%能量守恒
R_m=R(1:num);
R_p=R(num+1:2*num);
v_m_p=-u_Minus_Plus(:,t)*1i*w;
t_m_p=t_Minus_Plus(:,t);
p_m_p=real((v_m_p.'*conj(t_m_p))*ElemLen_m_power/2);
for iii=1:num
v_m(:,iii)=(-u_Minus(:,iii)*1i*w)*R_m(iii);
v_p(:,iii)=(-u_Plus(:,iii)*1i*w)*R_p(iii);
t_m(:,iii)=t_Minus(:,iii)*R_m(iii);
t_p(:,iii)=t_Plus(:,iii)*R_p(iii);
p_m(iii)=real((v_m(:,iii).'*conj(t_m(:,iii)))*ElemLen_m_power/2);
p_p(iii)=real((v_p(:,iii).'*conj(t_p(:,iii)))*ElemLen_p_power/2);
end











