function [MidPoint,midPoint,ElemLen,elemLen,NormalVector,normalVector,arrIndex,MarkB,...
    LenExtend,a,R,r,ElemLen_p_power]=ElemMesh_MainFrame_LS_Green_u_tot(H,WavLen,Nelem,Nwl,R_flaw,d_flaw,xC_flaw)
% 弯曲带缺陷
% 完好板构型% 中的格林函数计算，保留缺陷形状作后续处理。
% H: 半板厚
% R: 缺陷半径
% d: 缺陷深度
% xC_flaw: 缺陷中心位置

LenExtend = Nwl*WavLen; %半上边界长：从坐标原点距离左右截断位置的长度
L_elem = WavLen/Nelem; % 单元的长度
nElemExtend = round(LenExtend/L_elem); % 半上边界长单元数
R = 10.0;%外圈半径
r = 8.0;%内圈半径
a = 1*pi/2;%弯曲角度
nElemExtend1 = a*R/L_elem;%上弯曲边界长
nElemExtend2 = a*r/L_elem;%下弯曲边界长
transfer_matrix=[cos(a),sin(a);
    -sin(a),cos(a)];%坐标转换矩阵
transfer_matrix_flaw=[cos(a/2),sin(a/2);
    -sin(a/2),cos(a/2)];%缺陷坐标转换矩阵


nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end

%% 定义缺陷
% 现在设为1段缺陷
L_half_flaw=(R_flaw^2-(R_flaw-d_flaw)^2)^0.5; % 缺陷腐蚀的原边界一半长度
% ElemExtend1 = round((LenExtend-l)/L_elem); % 半下边界长单元数

L_flaw_range = abs(xC_flaw) + L_half_flaw;

MidPoint = [];
ElemLen = [];
NormalVector = [];
midPoint = [];
elemLen = [];
normalVector = [];
arrIndex = []; % 主框架四段，缺陷段，左右虚拟边界段 标志

MarkB = zeros(2,2); % 辅助边界在主框架中的位置
% 列代表每段边界


%%含弯曲节点板结构无缺陷边界
x1 = linspace(-LenExtend,0,nElemExtend+1); 
y1 = R*H*ones(1,nElemExtend+1); 
ListPointsTemp = [x1;y1].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N1 = max(size(MidPoint));
arrIndex = [arrIndex;N1];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];

t1_arc = linspace(0,a,nElemExtend1+1);
x1_arc = sin(t1_arc).*R;
y1_arc = cos(t1_arc).*R;
ListPointsTemp = [x1_arc;y1_arc].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N2 = max(size(MidPoint));
arrIndex = [arrIndex;N2];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];

x2 = linspace(0,LenExtend,nElemExtend+1); 
y2 = R*H*ones(1,nElemExtend+1);
coordinate1=[x2;
    y2];
ListPointsTemp = mtimes(transfer_matrix,coordinate1).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N3 = max(size(MidPoint));
arrIndex = [arrIndex;N3];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];


x3 = linspace(-LenExtend,0,nElemExtend+1); 
y3 = r*H*ones(1,nElemExtend+1); 
ListPointsTemp = [x3;y3].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N4 = max(size(MidPoint));
arrIndex = [arrIndex;N4];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];

x2_arc = sin(t1_arc).*r;
y2_arc = cos(t1_arc).*r;
ListPointsTemp = [x2_arc;y2_arc].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N5 = max(size(MidPoint));
arrIndex = [arrIndex;N5];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];

x4 = linspace(0,LenExtend,nElemExtend+1); 
y4 = r*H*ones(1,nElemExtend+1); 
coordinate2=[x4;
    y4];
ListPointsTemp = mtimes(transfer_matrix,coordinate2).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N6 = max(size(MidPoint));
arrIndex = [arrIndex;N6];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];


%% 缺陷
%{
L_halfArc = acos((R_flaw-d_flaw)/R_flaw)*R_flaw; % 裂纹弧长的一半

nElemHalfArc = round(L_halfArc/L_elem); % 缺陷上半部单元数的一半
if nElemHalfArc<5
    nElemHalfArc = 5;
end

t = acos((R_flaw-d_flaw)/R_flaw); % 半个弧度的角度
dt = t/nElemHalfArc;
t_arc = linspace(-t,t,2*nElemHalfArc+1);
x_arc = sin(t_arc).*R_flaw;
y_arc = ((cos(t_arc).*R_flaw+((r^2-(R_flaw^2-(R_flaw-d_flaw)^2))^0.5 -R_flaw+d_flaw)));
coordinate_flaw = [x_arc;
    y_arc];
ListPointsTemp = mtimes(transfer_matrix_flaw,coordinate_flaw).'; % 缺陷弧部 坐标
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NF = max(size(MidPoint));
arrIndex = [arrIndex;NF];
%}


% 因缺陷存在调整原边界与缺陷交接处单元，暂无效化。
%{
nL = Mark(1,1);
Point1 = MidPoint(nL,:) + 0.5*ElemLen(nL)*[-NormalVector(nL,2) NormalVector(nL,1)];
Point2 = [xC_flaw-L_half_flaw -H];
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo([Point1;Point2]);
MidPoint(nL,:) = MidPointTemp;
ElemLen(nL,:) = ElemLenTemp;
NormalVector(nL,:) = NormalVectorTemp;

nR = Mark(2,1);
Point1 = [xC_flaw+L_half_flaw -H];
Point2 = MidPoint(nR,:) + 0.5*ElemLen(nR)*[NormalVector(nR,2) -NormalVector(nR,1)];
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo([Point1;Point2]);
MidPoint(nR,:) = MidPointTemp;
ElemLen(nR,:) = ElemLenTemp;
NormalVector(nR,:) = NormalVectorTemp;
%}

MarkB(1,1) = arrIndex(1)-11;
MarkB(1,2) = arrIndex(2)+12;
MarkB(2,1) = arrIndex(4)-11;
MarkB(2,2) = arrIndex(5)+12;

nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end

xL1B = MidPoint(MarkB(1,1),1)+ElemLen(MarkB(1,1))/2;
xR1B = MidPoint(MarkB(1,2),1)*cos(-a)+MidPoint(MarkB(1,2),2)*sin(-a)-ElemLen(MarkB(1,2))/2;

%%左虚拟边界
x_L1 = ones(1,nElemBeneath_B+1)*xL1B;
y_L1 = linspace(r*H,R*H,nElemBeneath_B+1);
ListPointsTemp = [x_L1;y_L1].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NL1 = max(size(MidPoint));
arrIndex = [arrIndex;NL1];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];

%右虚拟边界
x_R1 = ones(1,nElemBeneath_B+1)*xR1B;
y_R1 = linspace(r*H,R*H,nElemBeneath_B+1);
coordinate_R1=[x_R1;
    y_R1];
ListPointsTemp = mtimes(transfer_matrix,coordinate_R1).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];

%%能量验证左虚拟边界
nElemBeneath_power=round(10*2*H/L_elem);
MarkB_(1,1) = 5;
MarkB_(1,2) = arrIndex(3)-5;
xL1B = MidPoint(MarkB_(1,1),1)+ElemLen(MarkB_(1,1))/2;
xR1B = MidPoint(MarkB_(1,2),1)*cos(-a)+MidPoint(MarkB_(1,2),2)*sin(-a)-ElemLen(MarkB_(1,2))/2;
x_L1 = ones(1,nElemBeneath_power+1)*xL1B;
y_L1 = linspace(r*H,R*H,nElemBeneath_power+1);
ListPointsTemp = [x_L1;y_L1].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NL1 = max(size(MidPoint));
arrIndex = [arrIndex;NL1];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];
%能量验证右虚拟边界
x_R1 = ones(1,nElemBeneath_power+1)*xR1B;
y_R1 = linspace(r*H,R*H,nElemBeneath_power+1);
ListPointsTemp = [x_R1;y_R1];
ListPointsTemp = mtimes(transfer_matrix,ListPointsTemp).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-a)+ListPointsTemp(:,2).*sin(-a);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-a)-ListPointsTemp(:,1).*sin(-a);
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
ElemLen_p_power=ElemLenTemp(1);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];
ListPointsTemp = ListPoints1;
[midPointTemp,elemLenTemp,normalVectorTemp]=MakeElemInfo(ListPointsTemp);
midPoint = [midPoint;midPointTemp];
elemLen = [elemLen;elemLenTemp];
normalVector = [normalVector;normalVectorTemp];
%{
plot(MidPoint(:,1),real(MidPoint(:,2)),'r')
hold on
%}