function [MidPoint,ElemLen,NormalVector]=MakeElemInfo(ListPoints)

MidPoint = [];
ElemLen = [];
NormalVector = [];

for t=1:max(size(ListPoints))-1
    midx1=(ListPoints(t+1,1)+ListPoints(t,1))/2;
    midy1=(ListPoints(t+1,2)+ListPoints(t,2))/2;
    CurMidPoint=[midx1, midy1];
    MidPoint=[MidPoint; CurMidPoint];
    L=((ListPoints(t+1,2)-ListPoints(t,2))^2+(ListPoints(t+1,1)-ListPoints(t,1))^2)^0.5;
    ElemLen=[ElemLen;L];
    normalvector1=-(ListPoints(t+1,2)-ListPoints(t,2))/L;
    normalvector2=+(ListPoints(t+1,1)-ListPoints(t,1))/L;
    CurNormalVector=[normalvector1, normalvector2];
    NormalVector=[NormalVector; CurNormalVector];
end