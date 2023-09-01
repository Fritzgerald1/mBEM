function [mat_r,mat_dr_D_dn,mat_r_,mat_dr_D_dn_,MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green_five(arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M 当前层的局部坐标
% 有缺陷
MidPoint_M = MidPoint(1:Ind(9),:);
ElemLen_M = ElemLen(1:Ind(9),:);
NormalVector_M = [NormalVector(1:Ind(3),:);-NormalVector(Ind(3)+1:Ind(7),:);-NormalVector(Ind(7)+1:Ind(8),:);NormalVector(Ind(8)+1:Ind(9),:)];  %法向：上表面向上，下表面向下
midPoint_M = midPoint(1:Ind(9),:);
elemLen_M = elemLen(1:Ind(9),:);
normalVector_M = [normalVector(1:Ind(3),:);-normalVector(Ind(3)+1:Ind(7),:);-normalVector(Ind(7)+1:Ind(8),:);normalVector(Ind(8)+1:Ind(9),:)];

Ind_M = [Ind(1:7,:); Ind(8);  Ind(9)];

[mat_r,mat_dr_D_dn] = GeneGeoInfoMat_MainFrame_five(arrX,MidPoint_M,ElemLen_M,NormalVector_M);
[mat_r_,mat_dr_D_dn_] = GeneGeoInfoMat_MainFrame_five(arrX,midPoint_M,elemLen_M,normalVector_M);