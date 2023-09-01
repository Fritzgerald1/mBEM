function [arrRes] = total_JS_Green4(MidPoint,NormalVector,Ind,MarkB,a,b,cl,ct,w,h1,rr, ...
    T,T_R,G,num_s,num,t,P_input,k_all)

%%%%%%%%%对应不同模态时，得到单位幅值Lamb波的位移和traction 
[nPt,temp]=size(MidPoint);

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
%%%%%在右半部分，各模态单位幅值Lamb波传播方向向右
x=MidPoint(:,1);
z=MidPoint(:,2)+rr;
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
ii=1:nPt;
u_Plus(2*ii-1,:)=u11(ii,:);
u_Plus(2*ii,:)=u22(ii,:);    
t_Plus(2*ii-1,:)=t11_tra(ii,:);
t_Plus(2*ii,:)=t22_tra(ii,:);
%% 修正项
Ind = 2*Ind;
MarkB = 2*MarkB;
nPtS = Ind(8); % nPtS: 表面单元数的二倍

matA_minus = zeros(nPtS,num);
matA_minus(1:MarkB(1,1),:) = -0.5*u_Minus(1:MarkB(1,1),:);
matA_minus(Ind(4)+1:MarkB(2,1),:) = -0.5*u_Minus(Ind(4)+1:MarkB(2,1),:);
matA_minus = matA_minus - T_R(1:nPtS,[1:MarkB(1,1) Ind(4)+1:MarkB(2,1) Ind(8)+1:Ind(9)])*u_Minus([1:MarkB(1,1) Ind(4)+1:MarkB(2,1) Ind(8)+1:Ind(9)],:);
matA_minus = matA_minus + G(1:nPtS,Ind(8)+1:Ind(9))*t_Minus(Ind(8)+1:Ind(9),:);

matA_plus = zeros(nPtS,num);
matA_plus(MarkB(1,2):Ind(4),:) = -0.5*u_Plus(MarkB(1,2):Ind(4),:);
matA_plus(MarkB(2,2):Ind(8),:) = -0.5*u_Plus(MarkB(2,2):Ind(8),:);
matA_plus = matA_plus - T_R(1:nPtS,[MarkB(1,2):Ind(4) MarkB(2,2):Ind(8) Ind(9)+1:Ind(10)])*u_Plus([MarkB(1,2):Ind(4) MarkB(2,2):Ind(8) Ind(9)+1:Ind(10)],:);
matA_plus = matA_plus + G(1:nPtS,Ind(9)+1:Ind(10))*t_Plus(Ind(9)+1:Ind(10),:);

matA = zeros(nPtS,1);
matA(1:MarkB(1,1),:) = -0.5*u_Minus_Plus(1:MarkB(1,1),t);
matA(Ind(4)+1:MarkB(2,1),:) = -0.5*u_Minus_Plus(Ind(4)+1:MarkB(2,1),t);
matA = matA - T_R(1:nPtS,[1:MarkB(1,1) Ind(4)+1:MarkB(2,1) Ind(8)+1:Ind(9)])*u_Minus_Plus([1:MarkB(1,1) Ind(4)+1:MarkB(2,1) Ind(8)+1:Ind(9)],t);
matA = matA + G(1:nPtS,Ind(8)+1:Ind(9))*t_Minus_Plus(Ind(8)+1:Ind(9),t);

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
matFinal(nPtS+n+1:nPtS+2*n,nPtS+n+1:nPtS+2*n) = u_plus(Ind(4)-n+1:Ind(4),:);
matFinal(nPtS+n+1:nPtS+2*n,Ind(4)-n+1:Ind(4)) = -eye(n);

%% 右端向量
matRight=zeros(nPtS+2*n,1);
matRight(1:nPtS,:) = - matA;
matRight(nPtS+1:nPtS+n,:) = - u_Minus_Plus(1:n,t);
arrRes=matFinal\matRight;%%%%大矩阵的运算
