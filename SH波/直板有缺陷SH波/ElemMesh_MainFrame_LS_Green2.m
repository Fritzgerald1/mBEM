function [MidPoint,ElemLen,NormalVector,arrIndex,L_flaw_range,L_half_flaw,Mark,MarkB,ElemLen_power]=ElemMesh_MainFrame_LS_Green2(H,xC_flaw,WavLen,Nelem,Nwl)

% ��ð幹���еĸ��ֺ������㣬����ȱ����״����������
% H: ����
% R: ȱ�ݰ뾶
% d: ȱ�����
% xC_flaw: ȱ������λ��

LenExtend = Nwl*WavLen; %���ϱ߽糤��������ԭ��������ҽض�λ�õĳ���
L_elem = WavLen/Nelem; % ��Ԫ�ĳ���
L_half_flaw=L_elem*22;
d = 0.5;
R=(d^2+L_half_flaw^2)/(2*d);
nElemExtend = round(LenExtend/L_elem); % ���ϱ߽糤��Ԫ��
nElemExtend_ = round((LenExtend-L_half_flaw)/L_elem); % ���±߽糤��Ԫ��
nElemBeneath_B=round(2*H/L_elem);
if nElemBeneath_B<5
    nElemBeneath_B = 5;
end

%% ����ȱ��
% ������Ϊ1��ȱ��


%L_half_flaw=(R^2-(R-d)^2)^0.5; % ȱ�ݸ�ʴ��ԭ�߽�һ�볤��
% ElemExtend1 = round((LenExtend-l)/L_elem); % ���±߽糤��Ԫ��

L_flaw_range = abs(xC_flaw) + L_half_flaw;

%%

MidPoint = [];
ElemLen = [];
NormalVector = [];
arrIndex = []; % ������ĶΣ�ȱ�ݶΣ���������߽�� ��־

% arrMark = zeros(4,1);

MarkB = zeros(2,2); % �����߽���������е�λ��
Mark = zeros(2,1); % ȱ�ݶ���������е�λ��
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
x3 = linspace(-LenExtend,-L_half_flaw,nElemExtend_+1); 
y3 = -H*ones(1,nElemExtend_+1); 
ListPointsTemp = [x3;y3].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N3 = max(size(MidPoint));
arrIndex = [arrIndex;N3];


x4 = linspace(L_half_flaw,LenExtend,nElemExtend_+1); 
y4 = -H*ones(1,nElemExtend_+1); 
ListPointsTemp = [x4;y4].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
N4 = max(size(MidPoint));
arrIndex = [arrIndex;N4];

%% 

MarkB(1,1) = find(MidPoint(1:arrIndex(1),1)<-L_flaw_range*2.0,1,'last');
MarkB(1,2) = arrIndex(1)+find(MidPoint(arrIndex(1)+1:arrIndex(2),1)>L_flaw_range*2.0,1,'first');
MarkB(2,1) = arrIndex(2)+find(MidPoint(arrIndex(2)+1:arrIndex(3),1)<-L_flaw_range*2.0,1,'last');
MarkB(2,2) = arrIndex(3)+find(MidPoint(arrIndex(3)+1:arrIndex(4),1)>L_flaw_range*2.0,1,'first');

%% ȱ��

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
%{
x5 = linspace(-L_half_flaw,L_half_flaw,50+1); 
y5 = -H*ones(1,50+1); 
ListPointsTemp = [x5;y5].';
[MidPointTemp,ElemLenTemp,NormalVectorTemp]=MakeElemInfo(ListPointsTemp);
MidPoint = [MidPoint;MidPointTemp];
ElemLen = [ElemLen;ElemLenTemp];
NormalVector = [NormalVector;NormalVectorTemp];
NF = max(size(MidPoint));
arrIndex = [arrIndex;NF];
%}
%%

% ��ǰȱ�����±���
Mark(1,1) = arrIndex(2)+find(MidPoint(arrIndex(2)+1:arrIndex(4),1) < xC_flaw-L_half_flaw,1,'last');
Mark(2,1) = arrIndex(3)+find(MidPoint(arrIndex(3)+1:arrIndex(4),1) > xC_flaw+L_half_flaw,1,'first');
%{
nL = Mark(1,1);
nR = Mark(2,1);
MidPoint = MidPoint([1:nL arrIndex(4)+1:arrIndex(5) nR:arrIndex(4)],:);
ElemLen = ElemLen([1:nL arrIndex(4)+1:arrIndex(5) nR:arrIndex(4)],:);
NormalVector = NormalVector([1:nL arrIndex(4)+1:arrIndex(5) nR:arrIndex(4)],:);
arrIndex = [arrIndex(1);arrIndex(2);nL+(arrIndex(5)-arrIndex(4))/2;arrIndex(5)-nR+nL+1];
%}
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
ElemLen_power = ElemLenTemp(1);
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
%{ 
plot(MidPoint(:,1),real(MidPoint(:,2)),'r')
hold on
%}