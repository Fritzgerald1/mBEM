%%%%%%%%%�ڽ�������ģ��ʱ�������ҽضϱ߽磬ȷ�����������߽�Ľ�λ��
clear
ff=[20e3,30e3,40e3,50e3,60e3,70e3,80e3,90e3,100e3,110e3,120e3,130e3,140e3,150e3,160e3,170e3,180e3,190e3,200e3,210e3,220e3,230e3,240e3,250e3,260e3,270e3,280e3,290e3,300e3];
%P_ref=zeros(29,2);
%P_tra=zeros(29,2);
Rtot = [];
nnn=1;
% for nnn=1:29
cl=6094;
ct=3263; 
f=ff(nnn);
w=2*pi*f;
rho=2700;
b=rho*ct^2;%nu����ģ��
a=rho*cl^2-2*b;%Lambda
E=b*(2*b+3*a)/(a+b);
v=a/(2*(a+b));
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
WavLen_S = 2*pi/k_all(t);

Nwl = 5;%%%%�����߽糡ȡ20������
Nelem = 100;%%%��Ԫ����Ϊ1/60������
L_elem = WavLen_A/Nelem;
LenExtend = Nwl*WavLen_A;
r = 50.8e-3;%��Ȧ�뾶
R = r+2*H;%��Ȧ�뾶
rr=0.5*(r+R);
alpha = 1*pi/2;%�����Ƕ�

theta = pi/4;%ȱ�ݽǶ�
deep = H/2;%ȱ�����
wide = 1*H/5;%ȱ�ݰ���

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
=ElemMesh_MainFrame_LS_Green_aocao(H,L_elem,LenExtend,r,R,alpha,theta,deep,wide);

[mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,...
MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green_aocao(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex);

[T,T_R,G]=GeneMat_aocao(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen_M,NormalVector_M,a,b,cl,ct,w);

[arrRes] = total_JS_Green_aocao(MidPoint_M,NormalVector_M,midPoint_M,normalVector_M,arrIndex_M,MarkB,a,b,cl,ct,w,H,rr,...
T,T_R,G,num_s,num,t,P_input,k_all,alpha);

%% 
% rrrr=arrIndex_M(8)+1:arrIndex_M(8)+2*num;
% R_=arrRes(2*rrrr-1);
rrrr=4*num;
R_=arrRes(end-rrrr:end);

% end
% [p_m,p_p,p_m_p]=power_1(a,b,cl,ct,w,H,rr,num_s,num,t,P_input,k_all,R_,ElemLen_m_power,...
% ElemLen_p_power,MidPoint_m_power,MidPoint_p_power);
% nnn=1;
% P_ref(nnn,:)=-p_m/p_m_p;
% P_tra(nnn,:)=p_p/p_m_p;
% disp(P_ref)
% disp(P_tra)
%disp(nnn);
%end
% [finish]=picture(MidPoint_M,arrRes,R_,arrIndex_M,P_ref,P_tra,R,r,num,ff);
% Finish=finish;