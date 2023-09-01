function [mat_r,mat_dr_D_dn]...
        =GeneGeoInfoMat_MainFrame_four(arrX,MidPoint,ElemLen,NormalVector)

nPoint = max(size(MidPoint));

mat_r = cell(nPoint,1);
mat_dr_D_dn = cell(nPoint,1);

AxisVec = [NormalVector(:,2), -NormalVector(:,1)]; %(n,2,1)
AxisVec = permute(AxisVec,[3,2,1]); %(1,2,n)
ElemLen = permute(ElemLen,[3,2,1]); %(n,1,1)->(1,1,n)
MidPoint = permute(MidPoint,[3,2,1]); %(n,2,1)->(1,2,n)
NormalVectorT = permute(NormalVector,[2,3,1]); %(n,2,1)->(2,1,n)


matCoordCur = bsxfun(@plus,MidPoint,bsxfun(@times,arrX.',bsxfun(@times,0.5*ElemLen,AxisVec)));

for loop1 = 1:nPoint %源点在各边界点中循环
    
    Coord_S = MidPoint(:,:,loop1); %源点坐标
    arrCoord_S = repmat(Coord_S,4,1);
    
    mat_dx_dy_temp = bsxfun(@minus,matCoordCur,arrCoord_S);
    mat_Cur_r = bsxfun(@hypot,mat_dx_dy_temp(:,1,:),mat_dx_dy_temp(:,2,:));
    mat_Cur_dr_D_dxi = mat_dx_dy_temp./repmat(mat_Cur_r,1,2,1);
    
    mat_Cur_dr_D_dn = bsxfun(@times,mat_Cur_dr_D_dxi(:,1,:),NormalVectorT(1,:,:)) + bsxfun(@times,mat_Cur_dr_D_dxi(:,2,:),NormalVectorT(2,:,:));
    
    mat_Cur_r = permute(mat_Cur_r,[1,3,2]);
    mat_Cur_dr_D_dn = permute(mat_Cur_dr_D_dn,[1,3,2]);
    
    mat_r(loop1) = {mat_Cur_r};
    mat_dr_D_dn(loop1) = {mat_Cur_dr_D_dn};
    
end



