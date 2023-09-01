%  freq10.0 每波长16个单元
% freq20.0 每波长8个单元


%定义介质参数

Rho = 1.0000;
Miu = 1.0000;

H = 1.0000; % 板厚2.0

curOmega = 20.0;

N_ind = 10;
arrN_ind = 0:1:2*N_ind-1;

cT = (Miu/Rho)^0.5;
kT = curOmega/cT;
arr_K1 = (kT^2 - (arrN_ind*pi/(2*H)).^2).^0.5;
arr_R1 = (arr_K1.^2 - kT^2).^0.5;

nPropMode = length(find(real(arr_K1)>0)); % 传播模态总数
N_ind_new = min(ceil(nPropMode/2+0.5),N_ind);
WavLen = 2*pi/arr_K1(1);


%%
xC_flaw = 0;
Nelem = 8; %每个波长的单元个数
Nwl = 20.125;

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



[MidPoint,ElemLen,NormalVector,arrIndex,L_flaw_range,L_half_flaw,Mark,MarkB,ElemLen_power]=ElemMesh_MainFrame_LS_Green2(H,xC_flaw,WavLen,Nelem,Nwl);

[mat_r,mat_dr_D_dn,...
MidPoint_M,ElemLen_M,NormalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green2(arrX,MidPoint,ElemLen,NormalVector,arrIndex);

[T,G]=GeneMat3(mat_r,mat_dr_D_dn,arrW,ElemLen_M,curOmega,Rho,Miu);

[arrRes,p,p_f_z] = total_JS_Green2(MidPoint_M,NormalVector_M,arrIndex_M,MarkB,Miu,Rho,H,curOmega,nPropMode,T,G,ElemLen_power);

x_plot_up = MidPoint(1:arrIndex_M(2),1);
x_plot_low = MidPoint([arrIndex_M(2)+1:arrIndex_M(3) arrIndex_M(4)+1:arrIndex_M(5) arrIndex_M(3)+1:arrIndex_M(4)],1);
nPropMode_R = 1:2*nPropMode;
nPropMode_I = 1:2*nPropMode;
figure(1)
plot(x_plot_up,real(arrRes(1:arrIndex_M(2),:)),'b')
hold on
figure(2)
plot(x_plot_up,imag(arrRes(1:arrIndex_M(2),:)),'b--')
hold 
figure(3)
plot(x_plot_low,real(arrRes([arrIndex_M(2)+1:arrIndex_M(3) arrIndex_M(4)+1:arrIndex_M(5) arrIndex_M(3)+1:arrIndex_M(4)],:)),'r')
hold on
figure(4)
plot(x_plot_low,imag(arrRes([arrIndex_M(2)+1:arrIndex_M(3) arrIndex_M(4)+1:arrIndex_M(5) arrIndex_M(3)+1:arrIndex_M(4)],:)),'r--')
hold 
figure(5)
plot(nPropMode_R,real(arrRes(arrIndex_M(5)+1:arrIndex_M(5)+2*nPropMode,:)),'r')
hold on
plot(nPropMode_I,imag(arrRes(arrIndex_M(5)+1:arrIndex_M(5)+2*nPropMode,:)),'r--')
hold 






