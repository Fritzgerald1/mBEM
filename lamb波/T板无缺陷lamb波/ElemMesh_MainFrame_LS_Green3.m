function [MidPoint,ElemLen,NormalVector,arrIndex,MarkB,...
    MidPoint_m_power,MidPoint_p_power,ElemLen_m_power,ElemLen_p_power,NormalVector_m_power,NormalVector_p_power]...
    =ElemMesh_MainFrame_LS_Green3(H,L_elem,LenExtend,r,alpha)
%%%%%%%%%%%%%边界元建模
nElemExtend = round(LenExtend/L_elem); % 半上边界长单元数
nElemExtend1 = round(alpha*r/L_elem);%上弯曲边界单元数
nElemExtend2 = round((LenExtend-r-H)/L_elem);

%%
MidPoint = [];
ElemLen = [];
NormalVector = [];
arrIndex = []; % 主框架四段，缺陷段，左右虚拟边界段 标志
MarkB = zeros(3,2); % 辅助边界在主框架中的位置
% 列代表每段边界

%%含弯曲节点板结构无缺陷边界
x1 = linspace(-LenExtend,0,nElemExtend+1); 
y1 = -H*ones(1,nElemExtend+1); 
ListPointsTemp = [x1;y1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N1 = max(size(MidPoint));
arrIndex = [arrIndex;N1];

x2 = linspace(0,LenExtend,nElemExtend+1); 
y2 = -H*ones(1,nElemExtend+1);
ListPointsTemp=[x2;y2].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N2 = max(size(MidPoint));
arrIndex = [arrIndex;N2];

x3 = linspace(-LenExtend,-r-H,nElemExtend2+1); 
y3 = H*ones(1,nElemExtend2+1); 
ListPointsTemp = [x3;y3].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N3 = max(size(MidPoint));
arrIndex = [arrIndex;N3];

t1_arc = linspace(2*alpha,alpha,nElemExtend1+1);
x1_arc = sin(t1_arc).*r-H-r;
y1_arc = cos(t1_arc).*r+H+r;
ListPointsTemp = [x1_arc;y1_arc].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N4 = max(size(MidPoint));
arrIndex = [arrIndex;N4];

y4 = linspace(r+H,LenExtend,nElemExtend2+1); 
x4 = -H*ones(1,nElemExtend2+1); 
ListPointsTemp = [x4;y4].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N5 = max(size(MidPoint));
arrIndex = [arrIndex;N5];

y5 = linspace(r+H,LenExtend,nElemExtend2+1); 
x5 = H*ones(1,nElemExtend2+1); 
ListPointsTemp = [x5;y5].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N6 = max(size(MidPoint));
arrIndex = [arrIndex;N6];

t2_arc = linspace(3*alpha,2*alpha,nElemExtend1+1);
x2_arc = sin(t2_arc).*r+r+H;
y2_arc = cos(t2_arc).*r+r+H;
ListPointsTemp = [x2_arc;y2_arc].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N7 = max(size(MidPoint));
arrIndex = [arrIndex;N7];

x6 = linspace(H+r,LenExtend,nElemExtend2+1); 
y6 = H*ones(1,nElemExtend2+1); 
ListPointsTemp=[x6;y6].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N8 = max(size(MidPoint));
arrIndex = [arrIndex;N8];

K=20;
KK=nElemExtend-nElemExtend2;
MarkB(1,1) = arrIndex(1)-KK-K;
MarkB(1,2) = arrIndex(1)+KK+K;
MarkB(2,1) = arrIndex(3)-K;
MarkB(2,2) = arrIndex(4)+K;
MarkB(3,1) = arrIndex(5)+K;
MarkB(3,2) = arrIndex(7)+K;

nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end

xL1B = MidPoint(MarkB(1,1),1)+ElemLen(MarkB(1,1))/2;
xR1B = MidPoint(MarkB(1,2),1)-ElemLen(MarkB(1,2))/2;
xR2B = MidPoint(MarkB(2,2),2)-ElemLen(MarkB(2,2))/2;

%%左虚拟边界
x_L1 = ones(1,nElemBeneath_B+1)*xL1B;
y_L1 = linspace(-H,H,nElemBeneath_B+1);
ListPointsTemp = [x_L1;y_L1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NL1 = max(size(MidPoint));
arrIndex = [arrIndex;NL1];

%右1虚拟边界
x_R1 = ones(1,nElemBeneath_B+1)*xR1B;
y_R1 = linspace(-H,H,nElemBeneath_B+1);
ListPointsTemp=[x_R1;y_R1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];

%右2虚拟边界
y_R2 = ones(1,nElemBeneath_B+1)*xR2B;
x_R2 = linspace(H,-H,nElemBeneath_B+1);
ListPointsTemp=[x_R2;y_R2].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];

%%能量验证左虚拟边界
nElemBeneath_power=round(2*H/L_elem);
if nElemBeneath_power<5
    nElemBeneath_power = 5;
end
MarkB_(1,1) = 1;
MarkB_(1,2) = arrIndex(2);
xL1B = MidPoint(MarkB_(1,1),1)-ElemLen(MarkB_(1,1))/2;
xR1B = MidPoint(MarkB_(1,2),1)*cos(-alpha)+MidPoint(MarkB_(1,2),2)*sin(-alpha)+ElemLen(MarkB_(1,2))/2;
x_L1 = ones(1,nElemBeneath_power+1)*xL1B;
y_L1 = linspace(-H,H,nElemBeneath_power+1);
ListPointsTemp = [x_L1;y_L1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
ElemLen_m_power=ElemLenTemp(1);           
MidPoint_m_power = MidPointTemp;
NormalVector_m_power = NormalVectorTemp;
%能量验证右虚拟边界
x_R1 = ones(1,nElemBeneath_power+1)*xR1B;
y_R1 = linspace(-H,H,nElemBeneath_power+1);
ListPointsTemp = [x_R1;y_R1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
ElemLen_p_power=ElemLenTemp(1);
MidPoint_p_power = MidPointTemp;
NormalVector_p_power = NormalVectorTemp;


figure(1)
plot(MidPoint(:,1),MidPoint(:,2),'b')
%}