function [finish]=picture_point(MidPoint,arrRes,Ind,R,r)
%%%%%%位移
u_outer=Ind(1)+1:Ind(2);
u1_outer=arrRes(2*u_outer-1);
u2_outer=arrRes(2*u_outer);
u_inner=Ind(4)+1:Ind(5);
u1_inner=arrRes(2*u_inner-1);
u2_inner=arrRes(2*u_inner);
xxxx1=asin(MidPoint(Ind(1)+1:Ind(2),1)/R);
xxxx2=asin(MidPoint(Ind(4)+1:Ind(5),1)/r);
figure(1)
plot(xxxx1,real(u1_outer),'b',xxxx1,imag(u1_outer),'b--',xxxx1,real(u2_outer),'r',xxxx1,imag(u2_outer),'r--')
ylabel('\it{u[outer wall]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{arcsin(x/R)}','FontName','Times New Roman','FontSize',20)
legend('u1-real','u1-imag','u2-real','u2-imag'); 
hold on
figure(2)
plot(xxxx2,real(u1_inner),'b',xxxx2,imag(u1_inner),'b--',xxxx2,real(u2_inner),'r',xxxx2,imag(u2_inner),'r--')
ylabel('\it{u[inner wall]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{arcsin(x/r)}','FontName','Times New Roman','FontSize',20)
legend('u1-real','u1-imag','u2-real','u2-imag'); 
hold on
finish=1;