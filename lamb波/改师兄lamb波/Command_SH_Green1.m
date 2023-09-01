% freq10.0 每波长16个单元
% freq20.0 每波长8个单元

%%%%%%%%%在建立整个模型时就不包含左右截断边界，确定了上下两边界的进位移，  
cl=5940;               
ct=3200; 
f=1e6;
w=2*pi*f;              
b=7800*ct^2;
a=7800*cl^2-2*b;
E=b*(2*b+3*a)/(a+b);
v=a/(2*(a+b));
H=0.5e-3;%%%%%%%%%%%二分之一的板厚
%{
%定义介质参数

Rho = 1.0000;
Miu = 1.0000;

H = 1.0000; % 板厚2.0

curOmega = 10.0;

N_ind = 10;
arrN_ind = 0:1:2*N_ind-1;

cT = (Miu/Rho)^0.5;
kT = curOmega/cT;
arr_K1 = (kT^2 - (arrN_ind*pi/(2*H)).^2).^0.5;
arr_R1 = (arr_K1.^2 - kT^2).^0.5;

nPropMode = length(find(real(arr_K1)>0)); % 传播模态总数
N_ind_new = min(ceil(nPropMode/2+0.5),N_ind);
WavLen = 2*pi/arr_K1(1);

Nelem = 16; %每个波长的单元个数
Nwl = 20.125;
%}
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


[MidPoint,ElemLen,NormalVector,arrIndex,MarkB,nElemExtend_h,nElemExtend_R]=ElemMesh_MainFrame_LS_Green_1(H);

[mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,...
MidPoint_M,ElemLen_M,NormalVector_M,arrIndex_M]...
=SUB_GeneGeoInfoNew_LS_Green1(arrX,MidPoint,ElemLen,NormalVector,arrIndex,nElemExtend_h,nElemExtend_R);

[T,T_R,G]=GeneMat3(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen_M,NormalVector_M,a,b,cl,ct,w);

[arrRes,num] = total_JS_Green1(MidPoint_M,NormalVector_M,arrIndex_M,MarkB,a,b,cl,ct,f,w,H,T,T_R,G);


x_plot_up = MidPoint(1+arrIndex_M(2):arrIndex_M(4),1);
x_plot_low = MidPoint(arrIndex_M(2)+1:arrIndex_M(4),1);
nPropMode_R = 1:4*num;
nPropMode_I = 1:4*num;
xx=arrIndex_M(2)+1:arrIndex_M(4);
figure(1)
plot(x_plot_up,real(arrRes(2*xx-1,:)),'b')
hold on
figure(2)
plot(x_plot_up,real(arrRes(2*xx,:)),'b')
hold on
%{
figure(2)
plot(x_plot_up,imag(arrRes(1:arrIndex_M(2),:)),'b--')
hold 
figure(3)
plot(x_plot_low,real(arrRes(arrIndex_M(2)+1:arrIndex_M(4),:)),'r')
hold on
figure(4)
plot(x_plot_low,imag(arrRes(arrIndex_M(2)+1:arrIndex_M(4),:)),'r--')
hold 

figure(5)
plot(nPropMode_R,real(arrRes(arrIndex_M(5)*2+1:arrIndex_M(5)*2+4*num,:)),'r')
hold on
plot(nPropMode_I,imag(arrRes(arrIndex_M(5)*2+1:arrIndex_M(5)*2+4*num,:)),'r--')
hold 
%}
