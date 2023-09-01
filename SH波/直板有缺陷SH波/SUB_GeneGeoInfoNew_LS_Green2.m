function [mat_r_up,mat_dr_D_dn_up,...
         MidPoint_M,ElemLen_M,NormalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green2...                                                                                     
         (arrX,MidPoint,ElemLen,NormalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M 当前层的局部坐标

%% 主框架

MidPoint_M = MidPoint(1:Ind(7),:);
ElemLen_M = ElemLen(1:Ind(7),:);
NormalVector_M = [NormalVector(1:Ind(2),:);-NormalVector(Ind(2)+1:Ind(5),:);NormalVector(Ind(5)+1:Ind(6),:);-NormalVector(Ind(6)+1:Ind(7),:)];  %法向：上表面向上，下表面向下

Ind_M = [Ind(1:5,:); Ind(6);  Ind(7)];

[mat_r_up,mat_dr_D_dn_up] = GeneGeoInfoMat_MainFrame_2(arrX,MidPoint_M,ElemLen_M,NormalVector_M);



