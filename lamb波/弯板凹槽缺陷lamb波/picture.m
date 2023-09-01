function[finish]=picture(MidPoint,arrRes,R_,arrIndex_M,P_ref,P_tra,R,r,num,ff)
%%%%%%%%%%位移U_tot
%{
u_outer=arrIndex_M(1)+1:arrIndex_M(2);
u_inner=[arrIndex_M(4)+1:arrIndex_M(5) arrIndex_M(7)+1:arrIndex_M(8) arrIndex_M(5)+1:arrIndex_M(6)];
u1_outer=abs(arrRes(2*u_outer-1));
u2_outer=abs(arrRes(2*u_outer));
u1_inner=abs(arrRes(2*u_inner-1));
u2_inner=abs(arrRes(2*u_inner));
xxxx1=asin(MidPoint(arrIndex_M(1)+1:arrIndex_M(2))/R);
xxxx2=asin(MidPoint([arrIndex_M(4)+1:arrIndex_M(5) arrIndex_M(7)+1:arrIndex_M(8) arrIndex_M(5)+1:arrIndex_M(6)])/r);
figure(1)
plot(xxxx1,real(u1_outer),'b',xxxx1,real(u2_outer),'r')
ylabel('\it{u(outer wall)}','FontName','Times New Roman','FontSize',20)
xlabel('\it{arcsin(x/R)}','FontName','Times New Roman','FontSize',20)
legend('u1','u2'); 
hold on
figure(2)
plot(xxxx2,real(u1_inner),'b',xxxx2,real(u2_inner),'r')
ylabel('\it{u(inner wall)}','FontName','Times New Roman','FontSize',20)
xlabel('\it{arcsin(x/r)}','FontName','Times New Roman','FontSize',20)
legend('u1','u2'); 
hold on
%}
%%%%%%%散射系数
nPropMode = 1:2*num;
R=abs(R_);
figure(3)
plot(nPropMode,real(R),'r')
hold on

%%%%%%%%反射、透射能量
%{
fff=ff(1:29)*1e-5;
figure(4)
plot(fff,real(P_ref(:,1)),'r')
axis([0.2 3 0 1])
ylabel('\it{S0 Mode Energy[R]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
figure(5)
plot(fff,real(P_ref(:,2)),'r')
axis([0.2 3 0 1])
ylabel('\it{A0 Mode Energy[R]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
figure(6)
plot(fff,real(P_tra(:,1)),'b')
axis([0.2 3 0 1])
ylabel('\it{S0 Mode Energy[T]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
figure(7)
plot(fff,real(P_tra(:,2)),'b')
axis([0.2 3 0 1])
ylabel('\it{A0 Mode Energy[T]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
%}
finish=1;
