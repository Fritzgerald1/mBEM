function [arrRes] = total_JS_Green3(MidPoint,NormalVector,midPoint,normalVector,Ind,MarkB,a,b,cl,ct,w,h1,rr, ...
    T,T_R,G,num_s,num,t,P_input,k_all,alpha)

%%%%%%%%%对应不同模态时，得到单位幅值Lamb波的位移和traction
[nPt,~]=size(MidPoint);

u_Minus=zeros(2*nPt,num);
u_Plus=zeros(2*nPt,num);
u_Minus_Plus=zeros(2*nPt,num);
t_Minus=zeros(2*nPt,num);
t_Plus=zeros(2*nPt,num);
t_Minus_Plus=zeros(2*nPt,num);
x=MidPoint(:,1);
z=MidPoint(:,2)-rr;
n1=NormalVector(:,1);
n2=NormalVector(:,2);
%%%%%在左半部分，各模态单位幅值Lamb波传播方向向左
for tt=1:num
    k_b=-k_all(tt);%%取负值
    p=(w^2/cl^2-k_b^2)^0.5;
    q=(w^2/ct^2-k_b^2)^0.5;
    if tt<=num_s
        P_xbs=Power_BS_ratio(w,cl,ct,a,b,-k_b,h1);
        P_xbs=double(P_xbs);
        A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比        
        u11(:,tt)=double((P_input/P_xbs)^0.5*(1i*k_b*cos(p*z)+q*A1*cos(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xbs)^0.5*(-p*sin(p*z)-1i*k_b*A1*sin(q*z)).*exp(1i*k_b*x));     
        Txx_tra=double((P_input/P_xbs)^0.5*(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z))).*exp(1i*k_b*x));
        Txz_tra=double((P_input/P_xbs)^0.5*(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z))).*exp(1i*k_b*x));
        Tzz_tra=double((P_input/P_xbs)^0.5*(-a*(k_b^2+p^2)*cos(p*z)-2*b*(p^2*cos(p*z)+1i*k_b*q*A1*cos(q*z))).*exp(1i*k_b*x));
        t11_tra(:,tt)=Txx_tra.*n1+Txz_tra.*n2;
        t22_tra(:,tt)=Txz_tra.*n1+Tzz_tra.*n2;
    elseif tt>num_s%%%%%反对称模态
        P_xba=Power_BA_ratio(w,cl,ct,a,b,-k_b,h1);
        P_xba=double(P_xba);
        A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比 
        u11(:,tt)=double((P_input/P_xba)^0.5*(1i*k_b*sin(p*z)-q*A2*sin(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xba)^0.5*(p*cos(p*z)-1i*k_b*A2*cos(q*z)).*exp(1i*k_b*x));
        Txx_tra=double((P_input/P_xba)^0.5*(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
        Txz_tra=double((P_input/P_xba)^0.5*(b*(2*1i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z))).*exp(1i*k_b*x));
        Tzz_tra=double((P_input/P_xba)^0.5*(-a*(k_b^2+p^2)*sin(p*z)-2*b*(p^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
        t11_tra(:,tt)=Txx_tra.*n1+Txz_tra.*n2;
        t22_tra(:,tt)=Txz_tra.*n1+Tzz_tra.*n2;
    end
end
ii=1:nPt;
u_Minus(2*ii-1,:)=u11(ii,:);
u_Minus(2*ii,:)=u22(ii,:);    
t_Minus(2*ii-1,:)=t11_tra(ii,:);
t_Minus(2*ii,:)=t22_tra(ii,:);
%%%%%在左半部分，入射模态Lamb波传播方向向右
for tt=1:num
    k_b=k_all(tt);%%取正值
    p=(w^2/cl^2-k_b^2)^0.5;
    q=(w^2/ct^2-k_b^2)^0.5;
    if tt<=num_s
        A1=2*1i*k_b*p*sin(p*h1)/((k_b^2-q^2)*sin(q*h1));%对称模式振幅比
        u11(:,tt)=(1i*k_b*cos(p*z)+q*A1*cos(q*z)).*exp(1i*k_b*x);
        u22(:,tt)=(-p*sin(p*z)-1i*k_b*A1*sin(q*z)).*exp(1i*k_b*x);
        Txx_tra=(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z))).*exp(1i*k_b*x);
        Txz_tra=(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z))).*exp(1i*k_b*x);
        Tzz_tra=(-a*(k_b^2+p^2)*cos(p*z)-2*b*(p^2*cos(p*z)+1i*k_b*q*A1*cos(q*z))).*exp(1i*k_b*x);
        t11_tra(:,tt)=Txx_tra.*n1+Txz_tra.*n2;
        t22_tra(:,tt)=Txz_tra.*n1+Tzz_tra.*n2;
    elseif tt>num_s%%%%%反对称模态
        A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比
        u11(:,tt)=(1i*k_b*sin(p*z)-q*A2*sin(q*z)).*exp(1i*k_b*x);
        u22(:,tt)=(p*cos(p*z)-1i*k_b*A2*cos(q*z)).*exp(1i*k_b*x);
        Txx_tra=(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x);
        Txz_tra=(b*(2*1i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z))).*exp(1i*k_b*x);
        Tzz_tra=(-a*(k_b^2+p^2)*sin(p*z)-2*b*(p^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x);
        t11_tra(:,tt)=Txx_tra.*n1+Txz_tra.*n2;
        t22_tra(:,tt)=Txz_tra.*n1+Tzz_tra.*n2;
    end
end
ii=1:nPt;
u_Minus_Plus(2*ii-1,:)=u11(ii,:);
u_Minus_Plus(2*ii,:)=u22(ii,:);    
t_Minus_Plus(2*ii-1,:)=t11_tra(ii,:);
t_Minus_Plus(2*ii,:)=t22_tra(ii,:);

x=midPoint(:,1);
z=midPoint(:,2)-rr;
n1=normalVector(:,1);
n2=normalVector(:,2);
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
        Txx_tra=double((P_input/P_xbs)^0.5*(-a*(p^2+k_b^2)*cos(p*z)+2*b*(-k_b^2*cos(p*z)+1i*k_b*A1*q*cos(q*z))).*exp(1i*k_b*x));
        Txz_tra=double((P_input/P_xbs)^0.5*(b*(-2*1i*k_b*p*sin(p*z)+(k_b^2-q^2)*A1*sin(q*z))).*exp(1i*k_b*x));
        Tzz_tra=double((P_input/P_xbs)^0.5*(-a*(k_b^2+p^2)*cos(p*z)-2*b*(p^2*cos(p*z)+1i*k_b*q*A1*cos(q*z))).*exp(1i*k_b*x));
        t11_tra(:,tt)=Txx_tra.*n1+Txz_tra.*n2;
        t22_tra(:,tt)=Txz_tra.*n1+Tzz_tra.*n2;
    elseif tt>num_s%%%%%反对称模态
        P_xba=Power_BA_ratio(w,cl,ct,a,b,k_b,h1);
        A2=-2*1i*k_b*p*cos(p*h1)/((k_b^2-q^2)*cos(q*h1));%反对称模式振幅比
        u11(:,tt)=double((P_input/P_xba)^0.5*(1i*k_b*sin(p*z)-q*A2*sin(q*z)).*exp(1i*k_b*x));
        u22(:,tt)=double((P_input/P_xba)^0.5*(p*cos(p*z)-1i*k_b*A2*cos(q*z)).*exp(1i*k_b*x));
        Txx_tra=double((P_input/P_xba)^0.5*(a*(-(k_b^2+p^2)*sin(p*z))+2*b*(-k_b^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
        Txz_tra=double((P_input/P_xba)^0.5*(b*(2*1i*k_b*p*cos(p*z)+(k_b^2-q^2)*A2*cos(q*z))).*exp(1i*k_b*x));
        Tzz_tra=double((P_input/P_xba)^0.5*(-a*(k_b^2+p^2)*sin(p*z)-2*b*(p^2*sin(p*z)-1i*k_b*q*A2*sin(q*z))).*exp(1i*k_b*x));
        t11_tra(:,tt)=Txx_tra.*n1+Txz_tra.*n2;
        t22_tra(:,tt)=Txz_tra.*n1+Tzz_tra.*n2;
    end
end
U11=u11*cos(alpha)+u22*sin(alpha);
U22=-u11*sin(alpha)+u22*cos(alpha);
T11_tra=t11_tra*cos(alpha)+t22_tra*sin(alpha);
T22_tra=-t11_tra*sin(alpha)+t22_tra*cos(alpha);
ii=1:nPt;
u_Plus(2*ii-1,:)=U11(ii,:);
u_Plus(2*ii,:)=U22(ii,:);    
t_Plus(2*ii-1,:)=T11_tra(ii,:);
t_Plus(2*ii,:)=T22_tra(ii,:);
%% 修正项
Ind = 2*Ind;
MarkB = 2*MarkB;
nPtS = Ind(6); % nPtS: 表面单元数的二倍

matA_minus = zeros(nPtS,num);
matA_minus(1:MarkB(1,1),:) = -0.5*u_Minus(1:MarkB(1,1),:);
matA_minus(Ind(3)+1:MarkB(2,1),:) = -0.5*u_Minus(Ind(3)+1:MarkB(2,1),:);
matA_minus = matA_minus - T_R(1:nPtS,[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)],:);
matA_minus = matA_minus + G(1:nPtS,Ind(6)+1:Ind(7))*t_Minus(Ind(6)+1:Ind(7),:);

matA_plus = zeros(nPtS,num);
matA_plus(MarkB(1,2):Ind(3),:) = -0.5*u_Plus(MarkB(1,2):Ind(3),:);
matA_plus(MarkB(2,2):Ind(6),:) = -0.5*u_Plus(MarkB(2,2):Ind(6),:);
matA_plus = matA_plus - T_R(1:nPtS,[MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)])*u_Plus([MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)],:);
matA_plus = matA_plus + G(1:nPtS,Ind(7)+1:Ind(8))*t_Plus(Ind(7)+1:Ind(8),:);

matA = zeros(nPtS,1);
u_Minus_Plus_tot = u_Minus_Plus(:,1)+u_Minus_Plus(:,2);
matA(1:MarkB(1,1),:) = -0.5*u_Minus_Plus_tot(1:MarkB(1,1));
matA(Ind(3)+1:MarkB(2,1),:) = -0.5*u_Minus_Plus_tot(Ind(3)+1:MarkB(2,1));
matA = matA - T_R(1:nPtS,[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus_Plus_tot([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)]);
matA = matA + G(1:nPtS,Ind(6)+1:Ind(7))*t_Minus_Plus(Ind(6)+1:Ind(7),t);

n=2*num;
matA_Minus=zeros(nPtS,n);
matA_Plus=zeros(nPtS,n);
u_minus=zeros(nPtS,n);
u_plus=zeros(nPtS,n);
tt=1:num;
ii=1:nPtS/2;
matA_Minus(2*ii-1,2*tt-1)=matA_minus(2*ii-1,tt);
matA_Minus(2*ii,2*tt)=matA_minus(2*ii,tt);
matA_Plus(2*ii-1,2*tt-1)=matA_plus(2*ii-1,tt);
matA_Plus(2*ii,2*tt)=matA_plus(2*ii,tt);
u_minus(2*ii-1,2*tt-1)=u_Minus(2*ii-1,tt);
u_minus(2*ii,2*tt)=u_Minus(2*ii,tt);    
u_plus(2*ii-1,2*tt-1)=u_Plus(2*ii-1,tt);
u_plus(2*ii,2*tt)=u_Plus(2*ii,tt);    

%% 大矩阵
matFinal = zeros(nPtS+2*n,nPtS+2*n);

% 修正项加入
matFinal(1:nPtS,1:nPtS) = T(1:nPtS,1:nPtS);
matFinal(1:nPtS,nPtS+1:nPtS+n) = matA_Minus;
matFinal(1:nPtS,nPtS+n+1:nPtS+2*n) = matA_Plus;

% 连续性条件加入
matFinal(nPtS+1:nPtS+n,nPtS+1:nPtS+n) = u_minus(1:n,:);
matFinal(nPtS+1:nPtS+n,1:n) = -eye(n);
matFinal(nPtS+n+1:nPtS+2*n,nPtS+n+1:nPtS+2*n) = u_plus(Ind(3)-n+1:Ind(3),:);
matFinal(nPtS+n+1:nPtS+2*n,Ind(3)-n+1:Ind(3)) = -eye(n);

%% 右端向量

matRight=zeros(nPtS+2*n,1);
matRight(1:nPtS,:) = - matA;
matRight(nPtS+1:nPtS+n,:) = - u_Minus_Plus_tot(1:n);
arrRes=matFinal\matRight;%%%%大矩阵的运算
%{
%%点源验证
matRight1=zeros(nPtS+2*n,nPtS+2*n);
matRight1(1:nPtS,1:nPtS) = G(1:nPtS,1:nPtS);
T_is = zeros(nPtS+2*n,1);
Ind=Ind/2;
T_is(2*(Ind(4)+((Ind(5)-Ind(4)+1)/2))-1,1)=1;
T_is(2*(Ind(4)+((Ind(5)-Ind(4)+1)/2)),1)=1;
matRight=matRight1*T_is;
arrRes=matFinal\matRight;%%%%大矩阵的运算
%}
