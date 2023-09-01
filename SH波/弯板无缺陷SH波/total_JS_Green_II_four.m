function [arrRes]...
         =total_JS_Green_II_four(MidPoint,NormalVector,midPoint,normalVector,Ind,cT,Miu,H,Omega,n,MarkB,T,G,rr)
    
% 现在的节点、单元长、法向全是以参与计算的当前层节点为准
% n: 传播模态数
X1=MidPoint(:,1);
X2=MidPoint(:,2)-rr;
X3=midPoint(:,1);
X4=midPoint(:,2)-rr;
[nPt,temp]=size(MidPoint);
       
u_Plus = zeros(nPt,n);
u_Minus = zeros(nPt,n);
u_minus_plus = zeros(nPt,n);
u_plus_minus = zeros(nPt,n);
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

        u_Stress_Plus_x1=Miu*1i*kesi*sin(beta*X4).*exp(1i*kesi*X3);
        u_Stress_Plus_x2=Miu*beta*cos(beta*X4).*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*normalVector(:,1)+u_Stress_Plus_x2.*normalVector(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*sin(beta*X2).*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=Miu*beta*cos(beta*X2).*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector(:,1)+u_Stress_Minua_x2.*NormalVector(:,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*sin(beta*X2).*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=Miu*beta*cos(beta*X2).*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1.*NormalVector(:,1)+u_Stress_minus_Plus_x2.*NormalVector(:,2);

    else
        u_Plus(:,loop4)=cos(beta*X4).*exp(1i*kesi*X3);
        u_Minus(:,loop4)=cos(beta*X2).*exp(-1i*kesi*X1);
        u_minus_plus(:,loop4)=cos(beta*X2).*exp(1i*kesi*X1);

        u_Stress_Plus_x1=Miu*1i*kesi*cos(beta*X4).*exp(1i*kesi*X3);
        u_Stress_Plus_x2=-Miu*beta*sin(beta*X4).*exp(1i*kesi*X3);
        t_Plus(:,loop4)=u_Stress_Plus_x1.*normalVector(:,1)+u_Stress_Plus_x2.*normalVector(:,2);

        u_Stress_Minua_x1=-Miu*1i*kesi*cos(beta*X2).*exp(-1i*kesi*X1);
        u_Stress_Minua_x2=-Miu*beta*sin(beta*X2).*exp(-1i*kesi*X1);
        t_Minus(:,loop4)= u_Stress_Minua_x1.*NormalVector(:,1)+u_Stress_Minua_x2.*NormalVector(:,2);

        u_Stress_minus_Plus_x1=Miu*1i*kesi*cos(beta*X2).*exp(1i*kesi*X1);
        u_Stress_minus_Plus_x2=-Miu*beta*sin(beta*X2).*exp(1i*kesi*X1);
        t_minus_plus(:,loop4)=u_Stress_minus_Plus_x1.*NormalVector(:,1)+u_Stress_minus_Plus_x2.*NormalVector(:,2);

    end

end

%% 修正项

nPtS = Ind(6); % nPtS: 表面单元数

matA_minus = zeros(nPtS,n);
matA_minus(1:MarkB(1,1),:) = -0.5*u_Minus(1:MarkB(1,1),:);
matA_minus(Ind(3)+1:MarkB(2,1),:) = -0.5*u_Minus(Ind(3)+1:MarkB(2,1),:);
matA_minus = matA_minus - T(1:nPtS,[1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)])*u_Minus([1:MarkB(1,1) Ind(3)+1:MarkB(2,1) Ind(6)+1:Ind(7)],:);
matA_minus = matA_minus + G(1:nPtS,Ind(6)+1:Ind(7))*t_Minus(Ind(6)+1:Ind(7),:);

matA_plus = zeros(nPtS,n);
matA_plus(MarkB(1,2):Ind(3),:) = -0.5*u_Plus(MarkB(1,2):Ind(3),:);
matA_plus(MarkB(2,2):Ind(6),:) = -0.5*u_Plus(MarkB(2,2):Ind(6),:);
matA_plus = matA_plus - T(1:nPtS,[MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)])*u_Plus([MarkB(1,2):Ind(3) MarkB(2,2):Ind(6) Ind(7)+1:Ind(8)],:);
matA_plus = matA_plus + G(1:nPtS,Ind(7)+1:Ind(8))*t_Plus(Ind(7)+1:Ind(8),:);

%% 大矩阵

matFinal = zeros(nPtS+2*n,nPtS+2*n);

% 修正项加入
matFinal(1:nPtS,1:nPtS) = T(1:nPtS,1:nPtS) + 0.5*eye(nPtS);
matFinal(1:nPtS,nPtS+1:nPtS+n) = matA_minus;
matFinal(1:nPtS,nPtS+n+1:nPtS+2*n) = matA_plus;

% 连续性条件加入
matFinal(nPtS+1:nPtS+n,nPtS+1:nPtS+n) = u_Minus(1:n,:);
matFinal(nPtS+1:nPtS+n,1:n) = -eye(n);
matFinal(nPtS+n+1:nPtS+2*n,nPtS+n+1:nPtS+2*n) = u_Plus(Ind(3)-n+1:Ind(3),:);
matFinal(nPtS+n+1:nPtS+2*n,Ind(3)-n+1:Ind(3)) = -eye(n);

%% 右端向量
Ap = 1;
%{
matA_minus_plus = zeros(nPtS,1);
matA_minus_plus(1:MarkB(1,1)) = -0.5*u_minus_plus(1:MarkB(1,1),1);
matA_minus_plus(Ind(3)+1:MarkB(2,1)) = -0.5*u_minus_plus(Ind(3)+1:MarkB(2,1),1);
matA_minus_plus = matA_minus_plus - T(1:nPtS,[1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)])*u_minus_plus([1:MarkB(1,1)  Ind(3)+1:MarkB(2,1)  Ind(6)+1:Ind(7)],1);
matA_minus_plus = matA_minus_plus + G(1:nPtS,Ind(6)+1:Ind(7))*t_minus_plus(Ind(6)+1:Ind(7),1);

matRight = zeros(nPtS+2*n,1);
matRight(1:nPtS) = - matA_minus_plus(:).*Ap;
matRight(nPtS+1:nPtS+n) = - u_minus_plus(1:n,1).*Ap;
%% 求解
arrRes=matFinal\matRight;
%}
%%点源验证
matRight1=zeros(nPtS+2*n,nPtS+2*n);
matRight1(1:nPtS,1:nPtS) = G(1:nPtS,1:nPtS);
T_is = zeros(nPtS+2*n,1);
T_is(Ind(4)+1,1)=1;
matRight=matRight1*T_is;
arrRes=matFinal\matRight;%%%%大矩阵的运算







