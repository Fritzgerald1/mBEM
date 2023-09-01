% freq10.0 每波长16个单元
% freq20.0 每波长8个单元


%定义介质参数

Rho = 1.0000;%密度
Miu = 1.0000;%泊松比

H = 1.0000; % 板厚2.0

curOmega = 10.0;

N_ind = 20;
arrN_ind = 0:1:2*N_ind-1;

cT = (Miu/Rho)^0.5;
kT = curOmega/cT;
arr_K1 = (kT^2 - (arrN_ind*pi/(2*H)).^2).^0.5;
arr_R1 = (arr_K1.^2 - kT^2).^0.5;

nPropMode = length(find(real(arr_K1)>0)); % 传播模态总数
N_ind_new = min(ceil(nPropMode/2+0.5),N_ind);
WavLen = 2*pi/arr_K1(1);

R_flaw = 1.0;
d_flaw = 0.5;
xC_flaw = 0;
Nelem = 16; %每个波长的单元个数
Nwl = 100.125;

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


[MidPoint,midPoint,ElemLen,elemLen,NormalVector,normalVector,arrIndex,MarkB,LenExtend,a,R,r,...
ElemLen_p_power]=ElemMesh_MainFrame_LS_Green_u_tot(H,WavLen,Nelem,Nwl,R_flaw,d_flaw,xC_flaw);

[mat_r,mat_dr_D_dn,mat_r_Dxi,mat_dr_Dxi,mat_r_,mat_dr_D_dn_,mat_r_Dxi_,mat_dr_Dxi_,MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green_four(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex);

[T,G,T_Dxi,G_Dxi]=GeneMat_four(mat_r,mat_dr_D_dn,mat_r_Dxi,mat_dr_Dxi,arrW,ElemLen_M,curOmega,Rho,Miu);

[t,g,t_Dxi,g_Dxi]=GeneMat_four(mat_r_,mat_dr_D_dn_,mat_r_Dxi_,mat_dr_Dxi_,arrW,elemLen_M,curOmega,Rho,Miu);

[arrRes]= total_JS_Green_II_four(MidPoint_M,NormalVector_M,midPoint_M,normalVector_M,...
arrIndex_M,Miu,Rho,H,curOmega,nPropMode,MarkB,T,G,t,g);

[U_tot_M,u_tot_M,U_tot_P,u_tot_P,T_tot_M,T_tot_P,t_tot_M,t_tot_P,U_ref_M,u_ref_M,T_ref_M,t_ref_M,u_inc_M,t_inc_M]=total_JS_Green_II_u_tot(MidPoint_M, ...
NormalVector_M,midPoint_M,normalVector_M,arrIndex_M,Miu,Rho,H,curOmega,nPropMode,MarkB,T,G,t,g,arrRes, ...
T_Dxi,G_Dxi,t_Dxi,g_Dxi);

%P功率检测
v_p=-U_tot_P*1i*curOmega;
p_p=(v_p.'*conj(T_tot_P)+conj(v_p.')*T_tot_P)*ElemLen_p_power/4;
vp=-u_tot_P*1i*curOmega;
pp=(vp.'*conj(t_tot_P)+conj(vp.')*t_tot_P)*ElemLen_p_power/4;

v_m=-U_ref_M*1i*curOmega
p_m=(v_m.'*conj(T_ref_M)+conj(v_m.')*T_ref_M)*ElemLen_p_power/4;
vm=-u_ref_M*1i*curOmega;
pm=(vm.'*conj(t_ref_M)+conj(vm.')*t_ref_M)*ElemLen_p_power/4;

x_plot_m = MidPoint(arrIndex(8)+1:arrIndex(9),2);
x_plot_p = midPoint(arrIndex(9)+1:arrIndex(10),2);

%U_tot
figure(1)
plot(x_plot_m,real(U_tot_M(:)),'b')
hold on
plot(x_plot_m,real(u_tot_M(:)),'r--')
hold on
figure(2)
plot(x_plot_m,imag(U_tot_M(:)),'b')
hold on
plot(x_plot_m,imag(u_tot_M(:)),'r--')
hold on
figure(3)
plot(x_plot_p,real(U_tot_P(:)),'b')
hold on
plot(x_plot_p,real(u_tot_P(:)),'r--')
hold on
figure(4)
plot(x_plot_p,imag(U_tot_P(:)),'b')
hold on
plot(x_plot_p,imag(u_tot_P(:)),'r--')
hold on


%T_tot
figure(5)
plot(x_plot_m,real(T_tot_M(:)),'b')
hold on
plot(x_plot_m,real(t_tot_M(:)),'r--')
hold on
figure(6)
plot(x_plot_m,imag(T_tot_M(:)),'b')
hold on
plot(x_plot_m,imag(t_tot_M(:)),'r--')
hold on

figure(7)
plot(x_plot_p,real(T_tot_P(:)),'b')
hold on
plot(x_plot_p,real(t_tot_P(:)),'r--')
hold on
figure(8)
plot(x_plot_p,imag(T_tot_P(:)),'b')
hold on
plot(x_plot_p,imag(t_tot_P(:)),'r--')
hold on
%}