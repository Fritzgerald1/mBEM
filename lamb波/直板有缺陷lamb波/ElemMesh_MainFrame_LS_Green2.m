function [MidPoint,ElemLen,NormalVector,arrIndex,MarkB]=ElemMesh_MainFrame_LS_Green2(H,L_elem,LenExtend,h,R)

% 完好板构型中的格林函数计算，保留缺陷形状作后续处理。
% H: 半板厚
% R: 缺陷半径
% d: 缺陷深度
% xC_flaw: 缺陷中心位置
%LenExtend = Nwl*WavLen; %半上边界长：从坐标原点距离左右截断位置的长度

nElemExtend_TOP = round(LenExtend/L_elem); % 半上边界长单元数
nElemExtend_BOTTOM = round((LenExtend-R)/L_elem);
%nElemExtend_h = round(h/L_elem1);
%nElemExtend_R = round(2*R/L_elem1);
nElemHalfArc = 5*round(pi*R/(2*L_elem));

nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end

%%

MidPoint = [];
ElemLen = [];
NormalVector = [];
arrIndex = []; % 主框架四段，缺陷段，左右虚拟边界段 标志
MarkB = zeros(2,2); % 辅助边界在主框架中的位置
% 列代表每段边界

%%
%上边界
x1 = linspace(-LenExtend,0,nElemExtend_TOP+1); 
y1 = H*ones(1,nElemExtend_TOP+1); 
ListPointsTemp = [x1;y1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N1 = max(size(MidPoint));
arrIndex = [arrIndex;N1];

x2 = linspace(0,LenExtend,nElemExtend_TOP+1); 
y2 = H*ones(1,nElemExtend_TOP+1); 
ListPointsTemp = [x2;y2].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N2 = max(size(MidPoint));
arrIndex = [arrIndex;N2];

%下边界
x3 = linspace(-LenExtend,-R,nElemExtend_BOTTOM+1); 
y3 = -H*ones(1,nElemExtend_BOTTOM+1); 
ListPointsTemp = [x3;y3].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N3 = max(size(MidPoint));
arrIndex = [arrIndex;N3];

x4 = linspace(R,LenExtend,nElemExtend_BOTTOM+1); 
y4 = -H*ones(1,nElemExtend_BOTTOM+1); 
ListPointsTemp = [x4;y4].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N4 = max(size(MidPoint));
arrIndex = [arrIndex;N4];

%%缺陷
%x5 = [-R*ones(1,nElemExtend_h+1),linspace(-R+L_elem1,R,nElemExtend_R),R*ones(1,nElemExtend_h)];
%y5 = [linspace(H,H-h,nElemExtend_h+1),(H-h)*ones(1,nElemExtend_R),linspace(H-h+L_elem1,H,nElemExtend_h)];
t = pi/2; % 半个弧度的角度
t_arc = linspace(-t,t,2*nElemHalfArc+1+1);
x5 = sin(t_arc).*R;
y5 = -(H-(cos(t_arc).*R));
ListPointsTemp = [x5;y5].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N5 = max(size(MidPoint));
arrIndex = [arrIndex;N5];
%{
%% 新增缺陷范围部左右竖向虚拟边界
MarkB(1,1) = arrIndex(1);
MarkB(2,1) = arrIndex(3)+15;
xL1B = MidPoint(MarkB(1,1),1)+ElemLen(MarkB(1,1))/2;

%左虚拟边界
x_L1 = ones(1,nElemBeneath_B-nElemExtend_h+1)*xL1B;
y_L1 = linspace(H-h,-H,nElemBeneath_B-nElemExtend_h+1);
ListPointsTemp = [x_L1;y_L1].';
%Inner_LR_Temp = [x_L1+L_elem; y_L1].';
%Inner_LR = [Inner_LR;Inner_LR_Temp];
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];

NL1 = max(size(MidPoint));
arrIndex = [arrIndex;NL1];

%右虚拟边界
x_R1 = ones(1,nElemBeneath_B-nElemExtend_h+1)*xL1B;
y_R1 = linspace(H-h,-H,nElemBeneath_B-nElemExtend_h+1);
ListPointsTemp = [x_R1;y_R1].';
%Inner_LR_Temp = [x_L1-L_elem; y_L1].';
%Inner_LR = [Inner_LR;Inner_LR_Temp];
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];

NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];
%}

MarkB(1,1) = arrIndex(1)-40;
MarkB(1,2) = arrIndex(1)+40;
MarkB(2,1) = arrIndex(3)-15;
MarkB(2,2) = arrIndex(4)+15;

%% 新增缺陷范围部左右竖向虚拟边界

xL1B = MidPoint(MarkB(1,1),1)+ElemLen(MarkB(1,1))/2;

xR1B = MidPoint(MarkB(1,2),1)-ElemLen(MarkB(1,2))/2;

%左虚拟边界
x_L1 = ones(1,nElemBeneath_B+1)*xL1B;
y_L1 = linspace(H,-H,nElemBeneath_B+1);
ListPointsTemp = [x_L1;y_L1].';
%Inner_LR_Temp = [x_L1+L_elem; y_L1].';
%Inner_LR = [Inner_LR;Inner_LR_Temp];
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];

NL1 = max(size(MidPoint));
arrIndex = [arrIndex;NL1];

%右虚拟边界
x_R1 = ones(1,nElemBeneath_B+1)*xR1B;
y_R1 = linspace(H,-H,nElemBeneath_B+1);
ListPointsTemp = [x_R1;y_R1].';
%Inner_LR_Temp = [x_L1-L_elem; y_L1].';
%Inner_LR = [Inner_LR;Inner_LR_Temp];
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];

NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];

%{
figure(1)
plot(MidPoint(:,1),real(MidPoint(:,2)),'b')
%}