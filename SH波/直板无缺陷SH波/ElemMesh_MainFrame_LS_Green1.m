function [MidPoint,ElemLen,NormalVector,arrIndex,MarkB,LenExtend,L_elem,ElemLen_p_power]=ElemMesh_MainFrame_LS_Green1(H,WavLen,Nelem,Nwl)

% ��ð幹���еĸ��ֺ������㣬����ȱ����״����������
% H: ����
% R: ȱ�ݰ뾶
% d: ȱ�����
% xC_flaw: ȱ������λ��

LenExtend = Nwl*WavLen; %���ϱ߽糤��������ԭ��������ҽض�λ�õĳ���
L_elem = WavLen/Nelem; % ��Ԫ�ĳ���
nElemExtend = round(LenExtend/L_elem); % ���ϱ߽糤��Ԫ��

nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end



%%

MidPoint = [];
ElemLen = [];
NormalVector = [];
arrIndex = []; % ������ĶΣ�ȱ�ݶΣ���������߽�� ��־
MarkB = zeros(2,2); % �����߽���������е�λ��
% �д���ÿ�α߽�

%%
%�ϱ߽�
x1 = linspace(-LenExtend,0,nElemExtend+1); 
y1 = H*ones(1,nElemExtend+1); 
ListPointsTemp = [x1;y1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N1 = max(size(MidPoint));
arrIndex = [arrIndex;N1];


x2 = linspace(0,LenExtend,nElemExtend+1); 
y2 = H*ones(1,nElemExtend+1); 
ListPointsTemp = [x2;y2].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N2 = max(size(MidPoint));
arrIndex = [arrIndex;N2];


%�±߽�
x3 = linspace(-LenExtend,0,nElemExtend+1); 
y3 = -H*ones(1,nElemExtend+1); 
ListPointsTemp = [x3;y3].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N3 = max(size(MidPoint));
arrIndex = [arrIndex;N3];


x4 = linspace(0,LenExtend,nElemExtend+1); 
y4 = -H*ones(1,nElemExtend+1); 
ListPointsTemp = [x4;y4].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N4 = max(size(MidPoint));
arrIndex = [arrIndex;N4];

%% 

MarkB(1,1) = 161-20;
MarkB(1,2) = 161+20;
MarkB(2,1) = 483-20;
MarkB(2,2) = 483+20;

%{ 
%% ����ȱ��
% ������Ϊ1��ȱ��
L_half_flaw=(R^2-(R-d)^2)^0.5; % ȱ�ݸ�ʴ��ԭ�߽�һ�볤��
% ElemExtend1 = round((LenExtend-l)/L_elem); % ���±߽糤��Ԫ��

L_flaw_range = abs(xC_flaw) + L_half_flaw;
ȱ��
L_halfArc = acos((R-d)/R)*R; % ���ƻ�����һ��

nElemHalfArc = round(L_halfArc/L_elem); % ȱ���ϰ벿��Ԫ����һ��
if nElemHalfArc<5
    nElemHalfArc = 5;
end

t = acos((R-d)/R); % ������ȵĽǶ�
dt = t/nElemHalfArc;
t_arc = linspace(-t,t,2*nElemHalfArc+1);
x_arc = sin(t_arc).*R;
y_arc = -(H-(cos(t_arc).*R-(R-d)));
ListPointsTemp = [x_arc;y_arc].'; % ȱ�ݻ��� ����
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NF = max(size(MidPoint));
arrIndex = [arrIndex;NF];
%}
%%
% ��ȱ�ݴ��ڵ���ԭ�߽���ȱ�ݽ��Ӵ���Ԫ������Ч����
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


%% ����ȱ�ݷ�Χ��������������߽�

xL1B = MidPoint(MarkB(1,1),1)+ElemLen(MarkB(1,1))/2;

xR1B = MidPoint(MarkB(1,2),1)-ElemLen(MarkB(1,2))/2;


%������߽�
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

%������߽�
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

%%������֤������߽�
nElemBeneath_power=round(10*2*H/L_elem);
MarkB_(1,1) = 5;
MarkB_(1,2) = arrIndex(2)-5;
xL1B = MidPoint(MarkB_(1,1),1)+ElemLen(MarkB_(1,1))/2;
xR1B = MidPoint(MarkB_(1,2),1)-ElemLen(MarkB_(1,2))/2;
x_L1 = ones(1,nElemBeneath_power+1)*xL1B;
y_L1 = linspace(-1*H,1*H,nElemBeneath_power+1);
ListPointsTemp = [x_L1;y_L1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NL1 = max(size(MidPoint));
arrIndex = [arrIndex;NL1];
%������֤������߽�
x_R1 = ones(1,nElemBeneath_power+1)*xR1B;
y_R1 = linspace(-1*H,1*H,nElemBeneath_power+1);
ListPointsTemp = [x_R1;y_R1].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
ElemLen_p_power=ElemLenTemp(1);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NR1 = max(size(MidPoint));
arrIndex = [arrIndex;NR1];
%{
plot(MidPoint(:,1),real(MidPoint(:,2)),'r')
hold on
%}