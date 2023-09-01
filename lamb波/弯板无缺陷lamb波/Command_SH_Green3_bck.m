%%%%%%%%%在建立整个模型时包含左右截断边界，确定了上下两边界的进位移
ff=[20e3,30e3,40e3,50e3,60e3,70e3,80e3,90e3,100e3,110e3,120e3,130e3,140e3,150e3,160e3,170e3,180e3,190e3,200e3,210e3,220e3,230e3,240e3,250e3,260e3,270e3,280e3,290e3,300e3];
%P_ref=zeros(29,2);
%P_tra=zeros(29,2);
%for nnn=1:29

%% 材料参数
cl=6094;               
ct=3263; 
f=300e3; % 频率
w=2*pi*f; % 角频率
rho=2700;
b=rho*ct^2; % 拉梅常数μ
a=rho*cl^2-2*b;% 拉梅常数λ
H=2.38e-3; % 二分之一的板厚
E=b*(2*b+3*a)/(a+b); % 杨氏模量
v=a/(2*(a+b)); % 泊松比

%% 频散曲线求波数
[k_s,k_a]=shiboshu(f);

%k_s=[0.374];
%k_a=[0.977];
%k_s=[0.760000];
%k_a=[1.61100];
num_s=length(k_s); % 对称模态数
k_all=[k_s,k_a]/H; % 波数
num=length(k_all); % S和A模态总数

t=2;%%%%% t确定入射模态,现在为A0模态入射
k=k_all(t);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 确定当前频率下的入射模态 ，现在入射A0模态

if t>num_s
    P_xa=Power_IA_ratio(w,cl,ct,a,b,k,H);%%%%%反对称模态入射波能流功率比
elseif t<=num_s
    P_xa=Power_IS_ratio(w,cl,ct,a,b,k,H);%%%%%对称模态入射波能流功率比
end
P_input=double(P_xa);
WavLen_A = 2*pi/k_all(num_s+1);
WavLen_S = 2*pi/k_all(1);

Nwl = 10;%%%%近场边界场取20个波长
Nelem = 60;%%%单元长度为1/60个波长
L_elem = WavLen_S/Nelem;
LenExtend = Nwl*WavLen_S;
r = 50.8e-3;%内圈半径
R = r+2*H;%外圈半径
rr=0.5*(r+R);
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

[MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex,MarkB,...
MidPoint_m_power,MidPoint_p_power,ElemLen_m_power,ElemLen_p_power,NormalVector_m_power,NormalVector_p_power]...
=ElemMesh_MainFrame_LS_Green3(H,L_elem,LenExtend,r,R,alpha);

[mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,...
MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green3(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex);

[T,T_R,G]=GeneMat3(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen_M,NormalVector_M,a,b,cl,ct,w);

[arrRes] = total_JS_Green3(MidPoint_M,NormalVector_M,midPoint_M,normalVector_M,arrIndex_M,MarkB,a,b,cl,ct,w,H,rr,...
T,T_R,G,num_s,num,t,P_input,k_all,alpha);
rrrr=arrIndex_M(6)+1:arrIndex_M(6)+2*num;
R_=arrRes(2*rrrr-1);
[finish]=picture_point(MidPoint_M,arrRes,arrIndex_M,R,r);
%[p_m,p_p,p_m_p]=power_1(a,b,cl,ct,w,H,rr,num_s,num,t,P_input,k_all,R_,ElemLen_m_power,...
%ElemLen_p_power,MidPoint_m_power,MidPoint_p_power);
%nnn=1;
%P_ref(nnn,:)=-p_m/p_m_p;
%P_tra(nnn,:)=p_p/p_m_p;
%disp(P_ref)
%disp(P_tra)
%disp(nnn);
%end
%[finish]=picture(MidPoint_M,arrRes,R_,arrIndex_M,P_ref,P_tra,R,r,num,ff);
Finish=finish;