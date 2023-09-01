function [mat_r,mat_dr_D_dn,MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green_four(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M 当前层的局部坐标
% 有缺陷
%{
MidPoint_M = [MidPoint([1:Ind(6) Ind(7)+1:Ind(9)],:)];
ElemLen_M = [ElemLen([1:Ind(6) Ind(7)+1:Ind(9)],:)];
NormalVector_M = [NormalVector(1:Ind(3),:);-NormalVector(Ind(3)+1:Ind(6),:);NormalVector(Ind(7)+1:Ind(8),:);-NormalVector(Ind(8)+1:Ind(9),:)];  %法向：上表面向上，下表面向下


Ind_M = [Ind(1:6,:); Ind(6)+Ind(8)-Ind(7);  Ind(6)+Ind(9)-Ind(7)];

[mat_r,mat_dr_D_dn] = GeneGeoInfoMat_MainFrame_four(arrX,MidPoint_M,ElemLen_M,NormalVector_M);
%}
%无缺陷
MidPoint_M = MidPoint(1:Ind(8),:);
ElemLen_M = ElemLen(1:Ind(8),:);
NormalVector_M = [NormalVector(1:Ind(3),:);-NormalVector(Ind(3)+1:Ind(6),:);-NormalVector(Ind(6)+1:Ind(7),:);NormalVector(Ind(7)+1:Ind(8),:)];  %法向：上表面向上，下表面向下
midPoint_M = midPoint(1:Ind(8),:);
elemLen_M = elemLen(1:Ind(8),:);
normalVector_M = [normalVector(1:Ind(3),:);-normalVector(Ind(3)+1:Ind(6),:);-normalVector(Ind(6)+1:Ind(7),:);normalVector(Ind(7)+1:Ind(8),:)];

Ind_M = [Ind(1:6,:); Ind(7);  Ind(8)];

[mat_r,mat_dr_D_dn] = GeneGeoInfoMat_MainFrame_four(arrX,MidPoint_M,ElemLen_M,NormalVector_M);
