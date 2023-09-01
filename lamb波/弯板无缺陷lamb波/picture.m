function[finish]=picture(MidPoint,arrRes,R_,arrIndex_M,P_ref,P_tra,R,r,num,ff)
%%%%%%位移
u_outer=arrIndex_M(1)+1:arrIndex_M(2);
u_inner=arrIndex_M(4)+1:arrIndex_M(5);
u1_outer=arrRes(2*u_outer-1);
u2_outer=arrRes(2*u_outer);
u1_inner=arrRes(2*u_inner-1);
u2_inner=arrRes(2*u_inner);
xxxx1=asin(MidPoint(arrIndex_M(1)+1:arrIndex_M(2))/R)/pi*2;
xxxx2=asin(MidPoint(arrIndex_M(4)+1:arrIndex_M(5))/r)/pi*2;

figure('Name','弯段外侧位移')
plot(xxxx1,real(u1_outer),'b',xxxx1,imag(u1_outer),'b--',xxxx1,real(u2_outer),'r',xxxx1,imag(u2_outer),'r--')
xlabel('\it{\theta/(\pi/2)}','FontName','Times New Roman','FontSize',20); ylabel('\it{u(outer wall)}','FontName','Times New Roman','FontSize',20);legend('u1-real','u1-imag','u2-real','u2-imag'); 
% xlabel('\it{arcsin(x/R)}','FontName','Times New Roman','FontSize',20)


figure('Name','弯段内侧位移')
plot(xxxx2,real(u1_inner),'b',xxxx2,imag(u1_inner),'b--',xxxx2,real(u2_inner),'r',xxxx2,imag(u2_inner),'r--')
xlabel('\it{\theta/(\pi/2)}','FontName','Times New Roman','FontSize',20); ylabel('\it{u(inner wall)}','FontName','Times New Roman','FontSize',20);legend('u1-real','u1-imag','u2-real','u2-imag'); legend('u1-real','u1-imag','u2-real','u2-imag'); 

%%%%%%%%散射系数
nPropMode = 1:2*num;
figure(3)
plot(nPropMode,real(R_(nPropMode)),'r.-')
hold on
plot(nPropMode,imag(R_(nPropMode)),'r--')
hold off
%%%%%%%%%能量
fff=ff(1:29)*1e-5;
figure(4)
plot(fff,real(P_ref(:,1)),'r')
% axis([0.2 3 0 1])
ylabel('\it{S0 Mode Energy[R]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
% hold off
figure(5)
plot(fff,real(P_ref(:,2)),'r')
% axis([0.2 3 0 1])
ylabel('\it{A0 Mode Energy[R]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
% hold off
figure(6)
plot(fff,real(P_tra(:,1)),'b')
% axis([0.2 3 0 1])
ylabel('\it{S0 Mode Energy[T]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
% hold off
figure(7)
plot(fff,real(P_tra(:,2)),'b')
% axis([0.2 3 0 1])
ylabel('\it{A0 Mode Energy[T]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
% hold off
finish=1;