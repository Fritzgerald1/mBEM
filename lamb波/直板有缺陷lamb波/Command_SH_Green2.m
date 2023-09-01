% freq10.0 ÿ����16����Ԫ
% freq20.0 ÿ����8����Ԫ

%%%%%%%%%�ڽ�������ģ��ʱ�Ͳ��������ҽضϱ߽磬ȷ�����������߽�Ľ�λ�ƣ�  
cl=1.8675;               
ct=1; 
rho=1;
f=0.1;
w=2*pi*f;              
b=rho*ct^2;
a=rho*cl^2-2*b;
H=1;%%%%%%%%%%%����֮һ�İ��

%[k_s,k_a]=shiboshu(f);
k_s=[0.374];
k_a=[0.977];
num_s=length(k_s);
k_all=[k_s,k_a]/H;
num=length(k_all);
t=2;%%%%%tȷ������ģ̬,����ΪA0ģ̬����
k=k_all(t);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȷ����ǰƵ���µ�����ģ̬ ����������A0ģ̬
if t>num_s
    P_xa=Power_IA_ratio(f,k,H);%%%%%���Գ�ģ̬���䲨�������ʱ�
elseif t<=num_s
    P_xa=Power_IS_ratio(f,k,H);%%%%%�Գ�ģ̬���䲨�������ʱ�
end
P_input=double(P_xa);
WavLen_A = 2*pi/k_all(num_s+1);
WavLen_S = 2*pi/k_all(1);
Nwl = 10;%%%%�����߽糡ȡ20������
Nelem = 12;%%%��Ԫ����Ϊ1/60������
L_elem = WavLen_A/Nelem;
LenExtend = Nwl*WavLen_A;
R=0.5;%%%%ȱ�ݶ���֮һ���
h=0.5;%%%%ȱ�����

nGauss = 4;
x1=-0.339981043584856;
x2=0.339981043584856;
x3=-0.861136311594053;
x4=0.861136311594053;
w1=0.652145154862546;
w2=0.652145154862546;
w3=0.347854845137454;
w4=0.347854845137454;
arrW = [w1,w2,w3,w4];
arrX = [x1,x2,x3,x4];

[MidPoint,ElemLen,NormalVector,arrIndex,MarkB]=ElemMesh_MainFrame_LS_Green2(H,L_elem,LenExtend,h,R);

[mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,...
MidPoint_M,ElemLen_M,NormalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green2(arrX,MidPoint,ElemLen,NormalVector,arrIndex);

[T,T_R,G]=GeneMat2(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen_M,NormalVector_M,a,b,cl,ct,w);

[arrRes] = total_JS_Green2(MidPoint_M,NormalVector_M,arrIndex_M,MarkB,a,b,cl,ct,f,w,H,T,T_R,G,num_s,num,t,P_input,k_all,k);

%[p_m,p_p,p_m_p]=power_1(arrRes,a,b,cl,ct,f,w,H,arrIndex_M,MidPoint);

x_plot_up = MidPoint(1:arrIndex_M(2),1);
x_plot_low = MidPoint([arrIndex_M(2)+1:arrIndex_M(3) arrIndex_M(4)+1:arrIndex_M(5) arrIndex_M(3)+1:arrIndex_M(4)],1);
nPropMode = 1:4*num;
tt=1:arrIndex_M(2);
u1=arrRes(2*tt-1);
u2=arrRes(2*tt);
rr=arrIndex_M(5)+1:arrIndex_M(5)+2*num;
R=arrRes(2*rr-1);
%{
figure(1)
plot(x_plot_up,real(arrRes(1:arrIndex_M(2),:)),'b')
hold on
figure(2)
plot(x_plot_up,imag(arrRes(1:arrIndex_M(2),:)),'b--')
hold 
figure(3)
plot(x_plot_low,real(arrRes(arrIndex_M(2)+1:arrIndex_M(4),:)),'r')
hold on
figure(4)
plot(x_plot_low,imag(arrRes(arrIndex_M(2)+1:arrIndex_M(4),:)),'r--')
%}
hold 
figure(1)
plot(x_plot_up,real(u1),'b')
hold on
figure(2)
plot(x_plot_up,real(u2),'r')
hold 

