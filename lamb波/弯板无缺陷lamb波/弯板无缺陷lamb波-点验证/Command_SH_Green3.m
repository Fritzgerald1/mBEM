%%%%%%%%%�ڽ�������ģ��ʱ�������ҽضϱ߽磬ȷ�����������߽�Ľ�λ��
ff=[110e3,120e3,130e3,140e3,150e3,160e3,170e3,180e3,190e3,200e3,210e3,220e3,230e3,240e3,250e3,260e3,270e3,280e3,290e3,300e3];
%P_ref=zeros(29,2);
%P_tra=zeros(29,2);
for nnn=1:2:length(ff)
cl=6094;               
ct=3263; 
f=200e3;
w=2*pi*f;
rho=2700;
b=rho*ct^2;%nu
a=rho*cl^2-2*b;%Lambda
H=2.38e-3;%%%%%%%%%%%����֮һ�İ��

[k_s,k_a]=shiboshu(f);
num_s=length(k_s);
k_all=[k_s,k_a]/H;
num=length(k_all);
t=2;%%%%%tȷ������ģ̬,����ΪA0ģ̬����
k=k_all(t);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȷ����ǰƵ���µ�����ģ̬ ����������A0ģ̬
if t>num_s
    P_xa=Power_IA_ratio(w,cl,ct,a,b,k,H);%%%%%���Գ�ģ̬���䲨�������ʱ�
elseif t<=num_s
    P_xa=Power_IS_ratio(w,cl,ct,a,b,k,H);%%%%%�Գ�ģ̬���䲨�������ʱ�
end
P_input=double(P_xa);
WavLen_A = 2*pi/k_all(num_s+1);
WavLen_S = 2*pi/k_all(1);
cp=w/k_all(1);
p=(w^2/cl^2-k^2)^0.5;
q=(w^2/ct^2-k^2)^0.5;
A1=2*1i*k*p*sin(p*H)/((k^2-q^2)*sin(q*H));
A2=-2*1i*k*p*cos(p*H)/((k^2-q^2)*cos(q*H));

Nwl = 10;%%%%�����߽糡ȡ20������
Nelem = 48;%%%��Ԫ����Ϊ1/60������
L_elem = WavLen_A/Nelem;
LenExtend = Nwl*WavLen_A;
r = 50.8e-3;%��Ȧ�뾶
R = r+2*H;%��Ȧ�뾶
rr=0.5*(r+R);
alpha = 1*pi/2;%�����Ƕ�

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

[MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex,MarkB,...
MidPoint_m_power,MidPoint_p_power,ElemLen_m_power,ElemLen_p_power,NormalVector_m_power,NormalVector_p_power]...
=ElemMesh_MainFrame_LS_Green3(H,L_elem,LenExtend,r,R,alpha);

[mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,...
MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green3(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex);

[T,T_R,G]=GeneMat3(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen_M,NormalVector_M,a,b,cl,ct,w);

[arrRes] = total_JS_Green3(MidPoint_M,NormalVector_M,midPoint_M,normalVector_M,arrIndex_M,MarkB,a,b,cl,ct,w,H,rr,...
T,T_R,G,num_s,num,t,P_input,k_all,alpha);
%% 
rrrr=2*arrIndex_M(6)+1;
R_=arrRes(rrrr:end);
Rtot = [Rtot, R_];
end
%% 
% 
% [finish]=picture_point(MidPoint_M,arrRes,arrIndex_M);
% 
% Finish=finish;