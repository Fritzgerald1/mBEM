 % freq10.0 每波长16个单元
% freq20.0 每波长8个单元


%定义介质参数

Rho = 1.0000;
Miu = 1.0000;

H = 1.0000; % 板厚2.0
%% 频散曲线
curOmega = 5; % 当前频率

N_ind = 10;
arrN_ind = 0:1:2*N_ind-1;
 
cT = (Miu/Rho)^0.5;
kT = curOmega/cT;
arr_K1 = (kT^2 - (arrN_ind*pi/(2*H)).^2).^0.5; % 波数
arr_R1 = (arr_K1.^2 - kT^2).^0.5;
%%
Prop_K1 = arr_K1(real(arr_K1)>0); % 传播模态波数
nPropMode = length(Prop_K1); % 传播模态总数
N_ind_new = min(ceil(nPropMode/2+0.5),N_ind);
WavLen = 2.*pi./arr_K1(1:2); % 当前模态波长


xC_flaw = 0;
Nelem = 8; %每个波长的单元个数
Nwl = 20.125; %每个单元的长度

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


[MidPoint,midPoint,ElemLen,elemLen,NormalVector,normalVector,arrIndex,MarkB,LenExtend,a,R,r,Mark,L_halfArc,...
ElemLen_m_power,ElemLen_p_power,MidPoint_m_power,MidPoint_p_power,NormalVector_m_power,NormalVector_p_power,...
ElemLen_power,scoop,ElemLen_scoop]=ElemMesh_MainFrame_LS_Green_five(H,WavLen,Nelem,Nwl,xC_flaw);

[mat_r,mat_dr_D_dn,mat_r_,mat_dr_D_dn_,MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green_five(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex);

[T,G]=GeneMat_five(mat_r,mat_dr_D_dn,arrW,ElemLen_M,curOmega,Rho,Miu);

[t,g]=GeneMat_five(mat_r_,mat_dr_D_dn_,arrW,elemLen_M,curOmega,Rho,Miu);
%% 

[arrRes,p,p_f_z]...
= total_JS_Green_II_five(MidPoint_M,NormalVector_M,midPoint_M,normalVector_M,arrIndex_M,Miu,Rho,H,curOmega,nPropMode,MarkB,T,G,t,g,Mark,ElemLen_power);

[P,p_m_p]=total_JS_Green_II_five_power(arrRes,arrIndex_M,Miu,Rho,H,curOmega,nPropMode,ElemLen_m_power, ...
ElemLen_p_power,MidPoint_m_power,MidPoint_p_power,NormalVector_m_power,NormalVector_p_power);

x_plot_up = linspace(-LenExtend-a*R/2,LenExtend+a*R/2,arrIndex_M(3));
x_plot_low = linspace(-LenExtend-a*r/2-L_halfArc+scoop*ElemLen_scoop/2,LenExtend+a*r/2+L_halfArc-scoop*ElemLen_scoop/2,arrIndex_M(7)-arrIndex_M(3));
nPropMode_R = 1:2*nPropMode;
nPropMode_I = 1:2*nPropMode;
figure(1)
plot(x_plot_up,real(arrRes(1:arrIndex_M(3),:)),'b')
hold on
figure(2)
plot(x_plot_up,imag(arrRes(1:arrIndex_M(3),:)),'b--')
hold on
figure(3)
plot(x_plot_low,real(arrRes([arrIndex_M(3)+1:Mark(1) arrIndex_M(6)+1:arrIndex_M(7) Mark(1)+1:arrIndex_M(6)],:)),'r')
hold on
figure(4)
plot(x_plot_low,imag(arrRes([arrIndex_M(3)+1:Mark(1) arrIndex_M(6)+1:arrIndex_M(7) Mark(1)+1:arrIndex_M(6)],:)),'r--')
hold on
figure(5)
plot(nPropMode_R,real(arrRes(arrIndex_M(7)+1:arrIndex_M(7)+2*nPropMode,:)),'r')
hold on
plot(nPropMode_I,imag(arrRes(arrIndex_M(7)+1:arrIndex_M(7)+2*nPropMode,:)),'r--')
hold on
