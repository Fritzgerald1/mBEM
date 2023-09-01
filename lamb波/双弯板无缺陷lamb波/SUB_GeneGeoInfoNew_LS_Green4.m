function [mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2,...
         MidPoint_M,ElemLen_M,NormalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green4...                                                                                     
         (arrX,MidPoint,ElemLen,NormalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M 当前层的局部坐标
MidPoint_M = MidPoint(1:Ind(10),:);
ElemLen_M = ElemLen(1:Ind(10),:);
NormalVector_M = [NormalVector(1:Ind(4),:);-NormalVector(Ind(4)+1:Ind(8),:);-NormalVector(Ind(8)+1:Ind(9),:);-NormalVector(Ind(9)+1:Ind(10),:)];  %法向：上表面向上，下表面向下

Ind_M = [Ind(1:8,:); Ind(9); Ind(10)];

[mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2] = GeneGeoInfoMat_MainFrame4(arrX,MidPoint_M,ElemLen_M,NormalVector_M);

