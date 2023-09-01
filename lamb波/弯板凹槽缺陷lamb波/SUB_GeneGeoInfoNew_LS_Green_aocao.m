function [mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2,...
         MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green_aocao...                                                                                     
         (arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M 当前层的局部坐标
MidPoint_M = MidPoint(1:Ind(10),:);
ElemLen_M = ElemLen(1:Ind(10),:);
NormalVector_M = [NormalVector(1:Ind(3),:);-NormalVector(Ind(3)+1:Ind(7),:);-NormalVector(Ind(7)+1:Ind(8),:);-NormalVector(Ind(8)+1:Ind(9),:);NormalVector(Ind(9)+1:Ind(10),:)];  %法向：上表面向上，下表面向下
midPoint_M = midPoint(1:Ind(10),:);
elemLen_M = elemLen(1:Ind(10),:);
normalVector_M = [normalVector(1:Ind(3),:);-normalVector(Ind(3)+1:Ind(7),:);-NormalVector(Ind(7)+1:Ind(8),:);-normalVector(Ind(8)+1:Ind(9),:);normalVector(Ind(9)+1:Ind(10),:)];  %法向：上表面向上，下表面向下

Ind_M = [Ind(1:8,:); Ind(9); Ind(10)];

[mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2] = GeneGeoInfoMat_MainFrame_aocao(arrX,MidPoint_M,ElemLen_M,NormalVector_M);


