function [finish]=picture_point(MidPoint,arrRes,Ind)
%%%%%%位移
u_outer=1:Ind(1);
u_inner=Ind(2)+1:Ind(3);
u1_outer=arrRes(2*u_outer-1)*1e3;
u2_outer=arrRes(2*u_outer)*1e3;
u1_inner=arrRes(2*u_inner-1)*1e3;
u2_inner=arrRes(2*u_inner)*1e3;
xxxx1=MidPoint(1:Ind(1),1)*1e3;
xxxx2=-MidPoint(Ind(2)+1:Ind(3),2)*1e3;
figure(1)
plot(xxxx1,real(u1_outer),'b')
ylabel('\it{u[mm]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{x[mm]}','FontName','Times New Roman','FontSize',20)
%legend('u1-real','u1-imag','u2-real','u2-imag'); 
hold on
figure(2)
plot(xxxx1,real(u2_outer),'r')
ylabel('\it{v[mm]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{y[mm]}','FontName','Times New Roman','FontSize',20)
%legend('u1-real','u1-imag','u2-real','u2-imag'); 
hold on
figure(3)
plot(xxxx2,real(u1_inner),'b')
ylabel('\it{u[mm]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{x[mm]}','FontName','Times New Roman','FontSize',20)
%legend('u1-real','u1-imag','u2-real','u2-imag'); 
hold on
figure(4)
plot(xxxx2,real(u2_inner),'r')
ylabel('\it{v[mm]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{y[mm]}','FontName','Times New Roman','FontSize',20)
%legend('u1-real','u1-imag','u2-real','u2-imag'); 
hold on

finish=1;