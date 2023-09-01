function [mat_r_up,mat_dr_D_dn_up,mat_r_Dxi,mat_dr_Dxi,...
         MidPoint_M,ElemLen_M,NormalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green1...                                                                                     
         (arrX,MidPoint,ElemLen,NormalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M ��ǰ��ľֲ�����

%% �����

MidPoint_M = MidPoint(1:Ind(8),:);
ElemLen_M = ElemLen(1:Ind(8),:);
NormalVector_M = [NormalVector(1:Ind(2),:);-NormalVector(Ind(2)+1:Ind(4),:);NormalVector(Ind(4)+1:Ind(5),:);-NormalVector(Ind(5)+1:Ind(6),:);-NormalVector(Ind(6)+1:Ind(7),:);NormalVector(Ind(7)+1:Ind(8),:)];  %�����ϱ������ϣ��±�������

Ind_M = [Ind(1:4,:); Ind(5); Ind(6); Ind(7); Ind(8)];

[mat_r_up,mat_dr_D_dn_up,mat_r_Dxi,mat_dr_Dxi] = GeneGeoInfoMat_MainFrame1(arrX,MidPoint_M,ElemLen_M,NormalVector_M);


