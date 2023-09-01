function [mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2,...
         MidPoint_M,ElemLen_M,NormalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green1...                                                                                     
         (arrX,MidPoint,ElemLen,NormalVector,Ind,nElemExtend_h,nElemExtend_R)

% MidPoint_M,ElemLen_M,NormalVector_M 当前层的局部坐标

%% 主框架

MidPoint_M = MidPoint(1:Ind(7),:);
ElemLen_M = ElemLen(1:Ind(7),:);
NormalVector_M = [NormalVector(1:Ind(2),:);-NormalVector(Ind(2)+1:Ind(4),:);NormalVector(Ind(4)+1:Ind(4)+nElemExtend_h,:);NormalVector(Ind(4)+1+nElemExtend_h:Ind(4)+nElemExtend_h+nElemExtend_R,:);NormalVector(Ind(4)+1+nElemExtend_h+nElemExtend_R:Ind(5),:);NormalVector(Ind(5)+1:Ind(6),:);-NormalVector(Ind(6)+1:Ind(7),:)];  %法向：上表面向上，下表面向下

Ind_M = [Ind(1:4,:); Ind(5); Ind(6); Ind(7)];

[mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2] = GeneGeoInfoMat_MainFrame1(arrX,MidPoint_M,ElemLen_M,NormalVector_M);


