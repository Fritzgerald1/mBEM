function [T,G]=GeneMat_five(mat_r,mat_dr_D_dn,arrW,ElemLen,Omega,Rho,Miu)

[nPoint,Temp] = size(ElemLen);


T = cell(nPoint,1);
G = cell(nPoint,1);


cT = (Miu/Rho)^0.5;
kT = Omega/cT;


for loop1 = 1:nPoint
    
    % 非对角元素
    % loop1
    
    cur_r = [mat_r{loop1}];
    cur_dr_D_dn = [mat_dr_D_dn{loop1}];
    
    MatT = -1i*kT/4*besselh(1,1,kT*cur_r).*cur_dr_D_dn;
    MatT(:,loop1)=0;
    MatT = arrW*MatT.*(0.5*ElemLen).';

    MatG = 1i/(4*Miu)*besselh(0,1,kT*cur_r);
    MatG(:,loop1)=0;
    MatG = arrW*MatG.*(0.5*ElemLen).';
    
    % 对角元素
    
    Lelem = ElemLen(loop1);
    GDiag_sta = -1/(2*pi*Miu)*Lelem*(log(Lelem/2)-1);
    %{
    U33_sta = -1/(2*pi*Miu)*log(r);
    U33 = 1i/(4*Miu)*besselh(0,1,kT*r);                         
    CurMatG = U33-U33_sta;
    %}
    matDiagElem_r = cur_r(:,loop1);
    MatGDiag_dyn = 1i/(4*Miu)*besselh(0,1,kT*matDiagElem_r) - (-1/(2*pi*Miu)*log(matDiagElem_r));
    GDiag_dyn = arrW*MatGDiag_dyn*(0.5*Lelem);
    MatG(:,loop1) = GDiag_sta + GDiag_dyn;
    
    T(loop1) = {MatT};
    G(loop1) = {MatG};

end

T = cell2mat(T);
G = cell2mat(G);

