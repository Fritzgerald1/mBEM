function [U_tot_M,u_tot_M,U_tot_P,u_tot_P,T_tot_M,T_tot_P,t_tot_M,t_tot_P,U_ref_M,u_ref_M,T_ref_M,t_ref_M,u_inc_M,t_inc_M]...
         =total_JS_Green_II_u_tot(MidPoint,NormalVector,midPoint,normalVector,Ind,Miu,Rho,H,Omega,n, ...
         MarkB,T,G,t,g,arrRes,T_Dxi,G_Dxi,t_Dxi,g_Dxi)
    
% 现在的节点、单元长、法向全是以参与计算的当前层节点为准
% n: 传播模态数
    
cT=(Miu/Rho)^0.5;

X1=MidPoint(:,1);
X2=MidPoint(:,2)-9;
X3=midPoint(:,1);
X4=midPoint(:,2)-9;
[nPt,temp]=size(MidPoint);
 
u_tot=arrRes(1:Ind(6));
R_M=arrRes(1+Ind(6):Ind(6)+n);
R_P=arrRes(1+n+Ind(6):Ind(6)+2*n);
       
u_Plus = zeros(nPt,n);
u_Minus = zeros(nPt,n);
u_minus_plus = zeros(nPt,n);
u_plus_minus = zeros(nPt,n);
u_Plus_dxi1=zeros(nPt,n);
u_Minus_dxi1=zeros(nPt,n);
u_minus_plus_dxi1=zeros(nPt,n);
u_plus_minus_dxi1=zeros(nPt,n);
u_Plus_dxi2=zeros(nPt,n);
u_Minus_dxi2=zeros(nPt,n);
u_minus_plus_dxi2=zeros(nPt,n);
u_plus_minus_dxi2=zeros(nPt,n);
t_Plus = zeros(nPt,n);
t_Minus = zeros(nPt,n);
t_minus_plus = zeros(nPt,n);
t_plus_minus = zeros(nPt,n);

for loop4=1:n
    nMode=loop4-1;
    xx=mod(nMode,2);
    beta=nMode*pi/(2*H);
    kesi=(Omega^2/cT^2-beta^2)^0.5;
    if xx==1
        u_Plus(:,loop4)=sin(beta*X4).*exp(1i*kesi*X3);
        u_Minus(:,loop4)=sin(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus(:,loop4)=sin(beta*X2).*exp(1i*kesi*X1);
        u_plus_minus(:,loop4)=sin(beta*X4).*exp(-1i*kesi*X3);
        u_Plus_dxi1(:,loop4)=1i*kesi*sin(beta*X4).*exp(1i*kesi*X3);
        u_Minus_dxi1(:,loop4)=-1i*kesi*sin(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus_dxi1(:,loop4)=1i*kesi*sin(beta*X2).*exp(1i*kesi*X1);
        u_plus_minus_dxi1(:,loop4)=-1i*kesi*sin(beta*X4).*exp(-1i*kesi*X3);
        u_Plus_dxi2(:,loop4)=beta*cos(beta*X4).*exp(1i*kesi*X3);
        u_Minus_dxi2(:,loop4)=beta*cos(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus_dxi2(:,loop4)=beta*cos(beta*X2).*exp(1i*kesi*X1);
        u_plus_minus_dxi2(:,loop4)=beta*cos(beta*X4).*exp(-1i*kesi*X2);

        u_Stress_Plus_x1=Miu*1i*kesi*sin(beta*X4).*exp(1i*kesi*X3);
        u_Stress_Plus_x2=Miu*beta*cos(beta*X4).*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*normalVector(:,1)+u_Stress_Plus_x2.*normalVector(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*sin(beta*X2).*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=Miu*beta*cos(beta*X2).*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector(:,1)+u_Stress_Minua_x2.*NormalVector(:,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*sin(beta*X2).*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=Miu*beta*cos(beta*X2).*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1.*NormalVector(:,1)+u_Stress_minus_Plus_x2.*NormalVector(:,2);

        u_Stress_Plus_minus_x1=-Miu*1i*kesi*sin(beta*X4).*exp(-1i*kesi*X3);
        u_Stress_Plus_minus_x2=Miu*beta*cos(beta*X4).*exp(-1i*kesi*X3);
        t_plus_minus(:,loop4)=u_Stress_Plus_minus_x1.*normalVector(:,1)+u_Stress_Plus_minus_x2.*normalVector(:,2);
    else
        u_Plus(:,loop4)=cos(beta*X4).*exp(1i*kesi*X3);
        u_Minus(:,loop4)=cos(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus(:,loop4)=cos(beta*X2).*exp(1i*kesi*X1);
        u_plus_minus(:,loop4)=cos(beta*X4).*exp(-1i*kesi*X3);
        u_Plus_dxi1(:,loop4)=1i*kesi*cos(beta*X4).*exp(1i*kesi*X3);
        u_Minus_dxi1(:,loop4)=-1i*kesi*cos(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus_dxi1(:,loop4)=1i*kesi*cos(beta*X2).*exp(1i*kesi*X1);
        u_plus_minus_dxi1(:,loop4)=-1i*kesi*cos(beta*X4).*exp(-1i*kesi*X3);
        u_Plus_dxi2(:,loop4)=-beta*sin(beta*X4).*exp(1i*kesi*X3);
        u_Minus_dxi2(:,loop4)=-beta*sin(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus_dxi2(:,loop4)=-beta*sin(beta*X2).*exp(1i*kesi*X1);
        u_plus_minus_dxi2(:,loop4)=-beta*sin(beta*X4).*exp(-1i*kesi*X3);

        u_Stress_Plus_x1=Miu*1i*kesi*cos(beta*X4).*exp(1i*kesi*X3);
        u_Stress_Plus_x2=-Miu*beta*sin(beta*X4).*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*normalVector(:,1)+u_Stress_Plus_x2.*normalVector(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*cos(beta*X2).*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=-Miu*beta*sin(beta*X2).*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector(:,1)+u_Stress_Minua_x2.*NormalVector(:,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*cos(beta*X2).*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=-Miu*beta*sin(beta*X2).*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1.*NormalVector(:,1)+u_Stress_minus_Plus_x2.*NormalVector(:,2);

        u_Stress_Plus_minus_x1=-Miu*1i*kesi*cos(beta*X4).*exp(-1i*kesi*X3);
        u_Stress_Plus_minus_x2=-Miu*beta*sin(beta*X4).*exp(-1i*kesi*X3);
        t_plus_minus(:,loop4)=u_Stress_Plus_minus_x1.*normalVector(:,1)+u_Stress_Plus_minus_x2.*normalVector(:,2);
    end

end
%uj对xi求导
u_Plus_dxi=u_Plus_dxi1.*normalVector(:,1)+u_Plus_dxi2.*normalVector(:,2);
u_Minus_dxi=u_Minus_dxi1.*NormalVector(:,1)+u_Minus_dxi2.*NormalVector(:,2);
u_minus_plus_dxi=u_minus_plus_dxi1.*NormalVector(:,1)+u_minus_plus_dxi2.*NormalVector(:,2);
u_plus_minus_dxi=u_plus_minus_dxi1.*normalVector(:,1)+u_plus_minus_dxi2.*normalVector(:,2);

%%计算全场位移
nPtS = Ind(9)-Ind(8); % nPtS: 表面单元数

matA_minus = zeros(nPtS,n);
matA_minus(:,:) = -u_Minus(Ind(8)+1:Ind(9),:);
matA_minus = matA_minus - T(Ind(8)+1:Ind(9),[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)],:);
matA_minus = matA_minus + G(Ind(8)+1:Ind(9),Ind(6)+1:Ind(7))*t_Minus(Ind(6)+1:Ind(7),:);

matA_plus = zeros(nPtS,n);
matA_plus = matA_plus - t(Ind(8)+1:Ind(9),[MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)])*u_Plus([MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)],:);
matA_plus = matA_plus + g(Ind(8)+1:Ind(9),Ind(7)+1:Ind(8))*t_Plus(Ind(7)+1:Ind(8),:);

matA_minus_plus = zeros(nPtS,1);
matA_minus_plus(:) = -u_minus_plus(Ind(8)+1:Ind(9),1);
matA_minus_plus = matA_minus_plus - T(Ind(8)+1:Ind(9),[1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)])*u_minus_plus([1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)],1);
matA_minus_plus = matA_minus_plus + G(Ind(8)+1:Ind(9),Ind(6)+1:Ind(7))*t_minus_plus(Ind(6)+1:Ind(7),1);

U_tot_M=-T(Ind(8)+1:Ind(9),1:Ind(6))*u_tot-matA_minus_plus-matA_minus*R_M-matA_plus*R_P;
U_ref_M=U_tot_M-u_minus_plus(Ind(8)+1:Ind(9),1);

nPtS = Ind(10)-Ind(9);

matA_minus = zeros(nPtS,n);
matA_minus = matA_minus - T(Ind(9)+1:Ind(10),[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)],:);
matA_minus = matA_minus + G(Ind(9)+1:Ind(10),Ind(6)+1:Ind(7))*t_Minus(Ind(6)+1:Ind(7),:);

matA_plus = zeros(nPtS,n);
matA_plus(:,:) = -u_Plus(Ind(9)+1:Ind(10),:);
matA_plus = matA_plus - t(Ind(9)+1:Ind(10),[MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)])*u_Plus([MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)],:);
matA_plus = matA_plus + g(Ind(9)+1:Ind(10),Ind(7)+1:Ind(8))*t_Plus(Ind(7)+1:Ind(8),:);

matA_minus_plus = zeros(nPtS,1);
matA_minus_plus = matA_minus_plus - T(Ind(9)+1:Ind(10),[1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)])*u_minus_plus([1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)],1);
matA_minus_plus = matA_minus_plus + G(Ind(9)+1:Ind(10),Ind(6)+1:Ind(7))*t_minus_plus(Ind(6)+1:Ind(7),1);

U_tot_P=-t(Ind(9)+1:Ind(10),1:Ind(6))*u_tot-matA_minus_plus-matA_minus*R_M-matA_plus*R_P;

u_tot_P=u_Plus(Ind(9)+1:Ind(10),:)*R_P;
u_ref_M=u_Minus(Ind(8)+1:Ind(9),:)*R_M;
u_inc_M=u_minus_plus(Ind(8)+1:Ind(9),1);
u_tot_M=u_ref_M+u_inc_M;

%%计算全场力
nPtS = Ind(9)-Ind(8);

matA_minus_dxi = zeros(nPtS,n);
matA_minus_dxi(:,:) = -u_Minus_dxi(Ind(8)+1:Ind(9),:);
matA_minus_dxi = matA_minus_dxi - T_Dxi(Ind(8)+1:Ind(9),[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)],:);
matA_minus_dxi = matA_minus_dxi + G_Dxi(Ind(8)+1:Ind(9),Ind(6)+1:Ind(7))*t_Minus(Ind(6)+1:Ind(7),:);

matA_plus_dxi = zeros(nPtS,n);
matA_plus_dxi = matA_plus_dxi - t_Dxi(Ind(8)+1:Ind(9),[MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)])*u_Plus([MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)],:);
matA_plus_dxi = matA_plus_dxi + g_Dxi(Ind(8)+1:Ind(9),Ind(7)+1:Ind(8))*t_Plus(Ind(7)+1:Ind(8),:);

matA_minus_plus_dxi = zeros(nPtS,1);
matA_minus_plus_dxi(:) = -u_minus_plus_dxi(Ind(8)+1:Ind(9),1);
matA_minus_plus_dxi = matA_minus_plus_dxi - T_Dxi(Ind(8)+1:Ind(9),[1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)])*u_minus_plus([1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)],1);
matA_minus_plus_dxi = matA_minus_plus_dxi + G_Dxi(Ind(8)+1:Ind(9),Ind(6)+1:Ind(7))*t_minus_plus(Ind(6)+1:Ind(7),1);

T_tot_M=-T_Dxi(Ind(8)+1:Ind(9),1:Ind(6))*u_tot-matA_minus_plus_dxi-matA_minus_dxi*R_M-matA_plus_dxi*R_P;
T_ref_M=T_tot_M-t_minus_plus(Ind(8)+1:Ind(9),1);

nPtS = Ind(10)-Ind(9);

matA_minus_dxi = zeros(nPtS,n);
matA_minus_dxi = matA_minus_dxi - T_Dxi(Ind(9)+1:Ind(10),[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)],:);
matA_minus_dxi = matA_minus_dxi + G_Dxi(Ind(9)+1:Ind(10),Ind(6)+1:Ind(7))*t_Minus(Ind(6)+1:Ind(7),:);

matA_plus_dxi = zeros(nPtS,n);
matA_plus_dxi(:,:) = -u_Plus_dxi(Ind(9)+1:Ind(10),:);
matA_plus_dxi = matA_plus_dxi - t_Dxi(Ind(9)+1:Ind(10),[MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)])*u_Plus([MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)],:);
matA_plus_dxi = matA_plus_dxi + g_Dxi(Ind(9)+1:Ind(10),Ind(7)+1:Ind(8))*t_Plus(Ind(7)+1:Ind(8),:);

matA_minus_plus_dxi = zeros(nPtS,1);
matA_minus_plus_dxi = matA_minus_plus_dxi - T_Dxi(Ind(9)+1:Ind(10),[1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)])*u_minus_plus([1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)],1);
matA_minus_plus_dxi = matA_minus_plus_dxi + G_Dxi(Ind(9)+1:Ind(10),Ind(6)+1:Ind(7))*t_minus_plus(Ind(6)+1:Ind(7),1);

T_tot_P=-T_Dxi(Ind(9)+1:Ind(10),1:Ind(6))*u_tot-matA_minus_plus_dxi-matA_minus_dxi*R_M-matA_plus_dxi*R_P;

t_tot_P=t_Plus(Ind(9)+1:Ind(10),:)*R_P;
t_ref_M=t_Minus(Ind(8)+1:Ind(9),:)*R_M;
t_inc_M=t_minus_plus(Ind(8)+1:Ind(9),1);
t_tot_M=t_ref_M+t_inc_M;
