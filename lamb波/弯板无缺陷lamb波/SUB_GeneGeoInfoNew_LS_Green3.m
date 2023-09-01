function [mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2,...
         MidPoint_M,ElemLen_M,NormalVector_M,midPoint_M,elemLen_M,normalVector_M,Ind_M]...
         =SUB_GeneGeoInfoNew_LS_Green3...                                                                                     
         (arrX,MidPoint,ElemLen,NormalVector,midPoint,elemLen,normalVector,Ind)

% MidPoint_M,ElemLen_M,NormalVector_M ��ǰ��ľֲ�����
MidPoint_M = MidPoint(1:Ind(8),:);
ElemLen_M = ElemLen(1:Ind(8),:);
NormalVector_M = [NormalVector(1:Ind(3),:);-NormalVector(Ind(3)+1:Ind(6),:);-NormalVector(Ind(6)+1:Ind(7),:);NormalVector(Ind(7)+1:Ind(8),:)];  %�����ϱ������ϣ��±������£����������ң��Ҹ���������
midPoint_M = midPoint(1:Ind(8),:);
elemLen_M = elemLen(1:Ind(8),:);
normalVector_M = [normalVector(1:Ind(3),:);-normalVector(Ind(3)+1:Ind(6),:);-normalVector(Ind(6)+1:Ind(7),:);normalVector(Ind(7)+1:Ind(8),:)];  %�����ϱ������ϣ��±������£����������ң��Ҹ���������

Ind_M = [Ind(1:6,:); Ind(7); Ind(8)];

[mat_r_up,mat_dr_D_dn_up,mat_r_x1,mat_r_x2] = GeneGeoInfoMat_MainFrame3(arrX,MidPoint_M,ElemLen_M,NormalVector_M);


