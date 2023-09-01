function [U_tot_M,u_tot_M,U_tot_P,u_tot_P,T_tot_M,T_tot_P,t_tot_M,t_tot_P,U_ref_M,u_ref_M,T_ref_M,t_ref_M,u_inc_M,t_inc_M]...
         =total_JS_Green1_u_tot(MidPoint,NormalVector,Ind,Miu,Rho,H,Omega,n, ...
         MarkB,T,G,arrRes,T_Dxi,G_Dxi)
    
% 现在的节点、单元长、法向全是以参与计算的当前层节点为准
% n: 传播模态数
    
cT=(Miu/Rho)^0.5;

X1=MidPoint(:,1);
X2=MidPoint(:,2);
[nPt,temp]=size(MidPoint);
 
u_tot=arrRes(1:Ind(4));
R_M=arrRes(1+Ind(4):Ind(4)+n);
R_P=arrRes(1+n+Ind(4):Ind(4)+2*n);
       
u_Plus = zeros(nPt,n);
u_Minus = zeros(nPt,n);
u_Plus_dxi1=zeros(nPt,n);
u_Minus_dxi1=zeros(nPt,n);
u_Plus_dxi2=zeros(nPt,n);
u_Minus_dxi2=zeros(nPt,n);
t_Plus = zeros(nPt,n);
t_Minus = zeros(nPt,n);

for loop4=1:n
    nMode=loop4-1;
    xx=mod(nMode,2);
    beta=nMode*pi/(2*H);
    kesi=(Omega^2/cT^2-beta^2)^0.5;
    if xx==1
        u_Plus(:,loop4)=sin(beta*X2).*exp(1i*kesi*X1);
        u_Minus(:,loop4)=sin(beta*X2).*exp(-1i*kesi*X1);
        u_Plus_dxi1(:,loop4)=1i*kesi*sin(beta*X2).*exp(1i*kesi*X1);
        u_Minus_dxi1(:,loop4)=-1i*kesi*sin(beta*X2).*exp(-1i*kesi*X1);
        u_Plus_dxi2(:,loop4)=beta*cos(beta*X2).*exp(1i*kesi*X1);
        u_Minus_dxi2(:,loop4)=beta*cos(beta*X2).*exp(-1i*kesi*X1);
        
        u_Stress_Plus_x1=Miu*1i*kesi*sin(beta*X2).*exp(1i*kesi*X1);
        u_Stress_Plus_x2=Miu*beta*cos(beta*X2).*exp(1i*kesi*X1);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*NormalVector(:,1)+u_Stress_Plus_x2.*NormalVector(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*sin(beta*X2).*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=Miu*beta*cos(beta*X2).*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector(:,1)+u_Stress_Minua_x2.*NormalVector(:,2);

    else
        u_Plus(:,loop4)=cos(beta*X2).*exp(1i*kesi*X1);
        u_Minus(:,loop4)=cos(beta*X2).*exp(-1i*kesi*X1);
        u_Plus_dxi1(:,loop4)=1i*kesi*cos(beta*X2).*exp(1i*kesi*X1);
        u_Minus_dxi1(:,loop4)=-1i*kesi*cos(beta*X2).*exp(-1i*kesi*X1);
        u_Plus_dxi2(:,loop4)=-beta*sin(beta*X2).*exp(1i*kesi*X1);
        u_Minus_dxi2(:,loop4)=-beta*sin(beta*X2).*exp(-1i*kesi*X1);

        u_Stress_Plus_x1=Miu*1i*kesi*cos(beta*X2).*exp(1i*kesi*X1);
        u_Stress_Plus_x2=-Miu*beta*sin(beta*X2).*exp(1i*kesi*X1);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*NormalVector(:,1)+u_Stress_Plus_x2.*NormalVector(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*cos(beta*X2).*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=-Miu*beta*sin(beta*X2).*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector(:,1)+u_Stress_Minua_x2.*NormalVector(:,2);

    end

end
%uj对xi求导
u_Plus_dxi=u_Plus_dxi1.*NormalVector(:,1)+u_Plus_dxi2.*NormalVector(:,2);
u_Minus_dxi=u_Minus_dxi1.*NormalVector(:,1)+u_Minus_dxi2.*NormalVector(:,2);

%%计算全场位移
nPtS = Ind(7)-Ind(6); % nPtS: 表面单元数

matA_minus = zeros(nPtS,n);
matA_minus(:,:) = -u_Minus(Ind(6)+1:Ind(7),:);
matA_minus = matA_minus - T(Ind(6)+1:Ind(7),[1:MarkB(1,1) Ind(2)+1:MarkB(2,1) Ind(4)+1:Ind(5)])*u_Minus([1:MarkB(1,1) Ind(2)+1:MarkB(2,1) Ind(4)+1:Ind(5)],:);
matA_minus = matA_minus + G(Ind(6)+1:Ind(7),Ind(4)+1:Ind(5))*t_Minus(Ind(4)+1:Ind(5),:);

matA_plus = zeros(nPtS,n);
matA_plus = matA_plus - T(Ind(6)+1:Ind(7),[MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)])*u_Plus([MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)],:);
matA_plus = matA_plus + G(Ind(6)+1:Ind(7),Ind(5)+1:Ind(6))*t_Plus(Ind(5)+1:Ind(6),:);

matA_minus_plus = zeros(nPtS,1);
matA_minus_plus(:) = -u_Plus(Ind(6)+1:Ind(7),1);
matA_minus_plus = matA_minus_plus - T(Ind(6)+1:Ind(7),[1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)])*u_Plus([1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)],1);
matA_minus_plus = matA_minus_plus + G(Ind(6)+1:Ind(7),Ind(4)+1:Ind(5))*t_Plus(Ind(4)+1:Ind(5),1);

U_tot_M=-T(Ind(6)+1:Ind(7),1:Ind(4))*u_tot-matA_minus_plus-matA_minus*R_M-matA_plus*R_P;
U_ref_M=U_tot_M-u_Plus(Ind(6)+1:Ind(7),1);

nPtS = Ind(8)-Ind(7);

matA_minus = zeros(nPtS,n);
matA_minus = matA_minus - T(Ind(7)+1:Ind(8),[1:MarkB(1,1) Ind(2)+1:MarkB(2,1) Ind(4)+1:Ind(5)])*u_Minus([1:MarkB(1,1) Ind(2)+1:MarkB(2,1) Ind(4)+1:Ind(5)],:);
matA_minus = matA_minus + G(Ind(7)+1:Ind(8),Ind(4)+1:Ind(5))*t_Minus(Ind(4)+1:Ind(5),:);

matA_plus = zeros(nPtS,n);
matA_plus(:,:) = -u_Plus(Ind(7)+1:Ind(8),:);
matA_plus = matA_plus - T(Ind(7)+1:Ind(8),[MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)])*u_Plus([MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)],:);
matA_plus = matA_plus + G(Ind(7)+1:Ind(8),Ind(5)+1:Ind(6))*t_Plus(Ind(5)+1:Ind(6),:);

matA_minus_plus = zeros(nPtS,1);
matA_minus_plus = matA_minus_plus - T(Ind(7)+1:Ind(8),[1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)])*u_Plus([1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)],1);
matA_minus_plus = matA_minus_plus + G(Ind(7)+1:Ind(8),Ind(4)+1:Ind(5))*t_Plus(Ind(4)+1:Ind(5),1);

U_tot_P=-T(Ind(7)+1:Ind(8),1:Ind(4))*u_tot-matA_minus_plus-matA_minus*R_M-matA_plus*R_P;

u_tot_P=u_Plus(Ind(7)+1:Ind(8),:)*R_P;
u_ref_M=u_Minus(Ind(6)+1:Ind(7),:)*R_M;
u_inc_M=u_Plus(Ind(6)+1:Ind(7),1);
u_tot_M=u_ref_M+u_inc_M;

%%计算全场力
nPtS = Ind(7)-Ind(6);

matA_minus_dxi = zeros(nPtS,n);
matA_minus_dxi(:,:) = -u_Minus_dxi(Ind(6)+1:Ind(7),:);
matA_minus_dxi = matA_minus_dxi - T_Dxi(Ind(6)+1:Ind(7),[1:MarkB(1,1) Ind(2)+1:MarkB(2,1) Ind(4)+1:Ind(5)])*u_Minus([1:MarkB(1,1) Ind(2)+1:MarkB(2,1) Ind(4)+1:Ind(5)],:);
matA_minus_dxi = matA_minus_dxi + G_Dxi(Ind(6)+1:Ind(7),Ind(4)+1:Ind(5))*t_Minus(Ind(4)+1:Ind(5),:);

matA_plus_dxi = zeros(nPtS,n);
matA_plus_dxi = matA_plus_dxi - T_Dxi(Ind(6)+1:Ind(7),[MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)])*u_Plus([MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)],:);
matA_plus_dxi = matA_plus_dxi + G_Dxi(Ind(6)+1:Ind(7),Ind(5)+1:Ind(6))*t_Plus(Ind(5)+1:Ind(6),:);

matA_minus_plus_dxi = zeros(nPtS,1);
matA_minus_plus_dxi(:) = -u_Plus_dxi(Ind(6)+1:Ind(7),1);
matA_minus_plus_dxi = matA_minus_plus_dxi - T_Dxi(Ind(6)+1:Ind(7),[1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)])*u_Plus([1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)],1);
matA_minus_plus_dxi = matA_minus_plus_dxi + G_Dxi(Ind(6)+1:Ind(7),Ind(4)+1:Ind(5))*t_Plus(Ind(4)+1:Ind(5),1);

T_tot_M=-T_Dxi(Ind(6)+1:Ind(7),1:Ind(4))*u_tot-matA_minus_plus_dxi-matA_minus_dxi*R_M-matA_plus_dxi*R_P;
T_ref_M=T_tot_M-t_Plus(Ind(6)+1:Ind(7),1);

nPtS = Ind(8)-Ind(7);

matA_minus_dxi = zeros(nPtS,n);
matA_minus_dxi = matA_minus_dxi - T_Dxi(Ind(7)+1:Ind(8),[1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)])*u_Minus([1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)],:);
matA_minus_dxi = matA_minus_dxi + G_Dxi(Ind(7)+1:Ind(8),Ind(4)+1:Ind(5))*t_Minus(Ind(4)+1:Ind(5),:);

matA_plus_dxi = zeros(nPtS,n);
matA_plus_dxi(:,:) = -u_Plus_dxi(Ind(7)+1:Ind(8),:);
matA_plus_dxi = matA_plus_dxi - T_Dxi(Ind(7)+1:Ind(8),[MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)])*u_Plus([MarkB(1,2):Ind(2) MarkB(2,2):Ind(4) Ind(5)+1:Ind(6)],:);
matA_plus_dxi = matA_plus_dxi + G_Dxi(Ind(7)+1:Ind(8),Ind(5)+1:Ind(6))*t_Plus(Ind(5)+1:Ind(6),:);

matA_minus_plus_dxi = zeros(nPtS,1);
matA_minus_plus_dxi = matA_minus_plus_dxi - T_Dxi(Ind(7)+1:Ind(8),[1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)])*u_Plus([1:MarkB(1,1)  Ind(2)+1:MarkB(2,1)  Ind(4)+1:Ind(5)],1);
matA_minus_plus_dxi = matA_minus_plus_dxi + G_Dxi(Ind(7)+1:Ind(8),Ind(4)+1:Ind(5))*t_Plus(Ind(4)+1:Ind(5),1);

T_tot_P=-T_Dxi(Ind(7)+1:Ind(8),1:Ind(4))*u_tot-matA_minus_plus_dxi-matA_minus_dxi*R_M-matA_plus_dxi*R_P;

t_tot_P=t_Plus(Ind(7)+1:Ind(8),:)*R_P;
t_ref_M=t_Minus(Ind(6)+1:Ind(7),:)*R_M;
t_inc_M=t_Plus(Ind(6)+1:Ind(7),1);
t_tot_M=t_ref_M+t_inc_M;