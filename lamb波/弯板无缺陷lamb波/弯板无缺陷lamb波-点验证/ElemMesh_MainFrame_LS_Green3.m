function [MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,arrIndex,MarkB,...
    MidPoint_m_power,MidPoint_p_power,ElemLen_m_power,ElemLen_p_power,NormalVector_m_power,NormalVector_p_power]...
    =ElemMesh_MainFrame_LS_Green3(H,L_elem,LenExtend,r,R,alpha)
%%%%%%%%%%%%%边界元建模
nElemExtend = round(LenExtend/L_elem); % 半上边界长单元数
nElemExtend1 = round(alpha*R/L_elem)-1;%上弯曲边界单元数
transfer_matrix=[cos(alpha),sin(alpha);
    -sin(alpha),cos(alpha)];%坐标转换矩阵

%%
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
y1 = R*ones(1,nElemExtend+1); 
ListPointsTemp = [x1;y1].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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

t1_arc = linspace(0,alpha,nElemExtend1+1);
x1_arc = sin(t1_arc).*R;
y1_arc = cos(t1_arc).*R;
ListPointsTemp = [x1_arc;y1_arc].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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
y2 = R*ones(1,nElemExtend+1);
coordinate1=[x2;
    y2];
ListPointsTemp = mtimes(transfer_matrix,coordinate1).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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
y3 = r*ones(1,nElemExtend+1); 
ListPointsTemp = [x3;y3].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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
y4 = r*ones(1,nElemExtend+1); 
coordinate2=[x4;
    y4];
ListPointsTemp = mtimes(transfer_matrix,coordinate2).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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

MarkB(1,1) = arrIndex(1)-11;
MarkB(1,2) = arrIndex(2)+12;
MarkB(2,1) = arrIndex(4)-11;
MarkB(2,2) = arrIndex(5)+12;

nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end

xL1B = MidPoint(MarkB(1,1),1)+ElemLen(MarkB(1,1))/2;
xR1B = MidPoint(MarkB(1,2),1)*cos(-alpha)+MidPoint(MarkB(1,2),2)*sin(-alpha)-ElemLen(MarkB(1,2))/2;

%%左虚拟边界
x_L1 = ones(1,nElemBeneath_B+1)*xL1B;
y_L1 = linspace(r,R,nElemBeneath_B+1);
ListPointsTemp = [x_L1;y_L1].';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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
y_R1 = linspace(r,R,nElemBeneath_B+1);
coordinate_R1=[x_R1;
    y_R1];
ListPointsTemp = mtimes(transfer_matrix,coordinate_R1).';
ListPoints1 = [];
ListPoints1(:,1)=ListPointsTemp(:,1).*cos(-alpha)+ListPointsTemp(:,2).*sin(-alpha);
ListPoints1(:,2)=ListPointsTemp(:,2).*cos(-alpha)-ListPointsTemp(:,1).*sin(-alpha);
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
nElemBeneath_power=round(2*H/L_elem);
if nElemBeneath_power<5
    nElemBeneath_power = 5;
end
MarkB_(1,1) = 1;
MarkB_(1,2) = arrIndex(3);
xL1B = MidPoint(MarkB_(1,1),1)-ElemLen(MarkB_(1,1))/2;
xR1B = MidPoint(MarkB_(1,2),1)*cos(-alpha)+MidPoint(MarkB_(1,2),2)*sin(-alpha)+ElemLen(MarkB_(1,2))/2;
x_L1 = ones(1,nElemBeneath_power+1)*xL1B;
y_L1 = linspace(r,R,nElemBeneath_power+1);
ListPointsTemp = [x_L1;y_L1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
ElemLen_m_power=ElemLenTemp(1);           
MidPoint_m_power = MidPointTemp;
NormalVector_m_power = NormalVectorTemp;
%能量验证右虚拟边界
x_R1 = ones(1,nElemBeneath_power+1)*xR1B;
y_R1 = linspace(r,R,nElemBeneath_power+1);
ListPointsTemp = [x_R1;y_R1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
ElemLen_p_power=ElemLenTemp(1);
MidPoint_p_power = MidPointTemp;
NormalVector_p_power = NormalVectorTemp;

%{
figure(1)
plot(MidPoint(:,1),real(MidPoint(:,2)),'b')
%}