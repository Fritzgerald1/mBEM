%%%%%%%%%在建立整个模型时包含左右截断边界，确定了上下两边界的进位移
ff=[100e3,110e3,120e3,130e3,140e3,150e3,160e3,170e3,180e3,190e3,200e3,210e3,220e3,230e3,240e3,250e3,260e3,270e3,280e3,290e3,300e3];
%P_ref=zeros(21,2);
%P_tra=zeros(21,2);
%for nnn=1:21
cl=6094;               
ct=3263; 
f=ff(1);
w=2*pi*f;
rho=2700;
b=rho*ct^2;%nu
a=rho*cl^2-2*b;%Lambda
H=1.5e-3;%%%%%%%%%%%二分之一的板厚

[k_s,k_a]=shiboshu(f);
num_s=length(k_s);
k_all=[k_s,k_a]/H;
num=length(k_all);
t=1;%%%%%t确定入射模态,现在为S0模态入射
k=k_all(t);%%%%确定当前频率下的入射模态 ，现在入射S0模态
if t>num_s
    P_xa=Power_IA_ratio(w,cl,ct,a,b,k,H);%%%%%反对称模态入射波能流功率比
elseif t<=num_s
    P_xa=Power_IS_ratio(w,cl,ct,a,b,k,H);%%%%%对称模态入射波能流功率比
end
P_input=double(P_xa);
WavLen = 2*pi/k_all(num_s+1);
WavLen_s= 2*pi/k_all(1);

Nwl = 10;%%%%近场边界场取20个波长
Nelem = 24;%%%单元长度为1/24个波长
L_elem = WavLen/Nelem;
LenExtend = Nwl*WavLen_s;
r = 50.8e-3;%内圈半径
R_a = r+2*H+3e-3;%外圈半径
R_b = r+2*H;
rr=0.5*(r+R_b);
alpha = 1*pi/2;%弯曲角度

%%%%nGauss = 4;
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

[MidPoint,ElemLen,NormalVector,arrIndex,MarkB,...
MidPoint_m_power,MidPoint_p_power,ElemLen_m_power,ElemLen_p_power,NormalVector_m_power,NormalVector_p_power]...
=ElemMesh_MainFrame_LS_Green4(H,L_elem,LenExtend,r,R_a,R_b,alpha);

[mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,...
MidPoint_M,ElemLen_M,NormalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green4(arrX,MidPoint,ElemLen,NormalVector,arrIndex);

[T,T_R,G]=GeneMat4(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen_M,NormalVector_M,a,b,cl,ct,w);

[arrRes] = total_JS_Green4(MidPoint_M,NormalVector_M,arrIndex_M,MarkB,a,b,cl,ct,w,H,rr,...
T,T_R,G,num_s,num,t,P_input,k_all);

rrrr=arrIndex_M(8)+1:arrIndex_M(8)+2*num;
R_=arrRes(2*rrrr-1);

[p_m,p_p,p_m_p]=power_1(a,b,cl,ct,w,H,rr,num_s,num,t,P_input,k_all,R_,ElemLen_m_power,...
ElemLen_p_power,MidPoint_m_power,MidPoint_p_power);
nnn=1;
P_ref(nnn,:)=-p_m/p_m_p;
P_tra(nnn,:)=p_p/p_m_p;
disp(nnn);
%end
%{
disp(P_ref)
disp(P_tra)
nPropMode = 1:2*num;
figure(5)
plot(nPropMode,real(R_),'r')
hold on
plot(nPropMode,imag(R_),'r--')
hold on
%}
fff=ff(1:21)*1e-5;
figure(1)
plot(fff,real(P_ref(:,1)),'r')
axis([1 3 0 1])
ylabel('\it{S0 Mode Energy[R]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
figure(2)
plot(fff,real(P_ref(:,2)),'r')
axis([1 3 0 1])
ylabel('\it{A0 Mode Energy[R]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
figure(3)
plot(fff,real(P_tra(:,1)),'b')
axis([1 3 0 1])
ylabel('\it{S0 Mode Energy[T]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
figure(4)
plot(fff,real(P_tra(:,2)),'b')
axis([1 3 0 1])
ylabel('\it{A0 Mode Energy[T]}','FontName','Times New Roman','FontSize',20)
xlabel('\it{f[100kHz]}','FontName','Times New Roman','FontSize',20)
set(gca,'FontSize',20)
hold on
