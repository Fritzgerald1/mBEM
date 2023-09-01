function [T,T_R,G]=GeneMat3(mat_r,mat_dr_D_dn,mat_r_x1,mat_r_x2,arrW,arrX,ElemLen,NormalVector,a,b,cl,ct,w)

k1=w/cl;
k2=w/ct;
oula=0.577215667901532860;
n=10;%%%%%%%取展开式的前10项
py_n=0;
for x=1:n
    py_n=py_n+1/x;
end
py_n1=py_n+1/(n+1);
[nPoint,Temp] = size(ElemLen);
n1=zeros(4,nPoint);
n2=zeros(4,nPoint);
for ii = 1:4
    n1(ii,:)=NormalVector(:,1);
    n2(ii,:)=NormalVector(:,2);
end
t11 = cell(nPoint,1);
t12 = cell(nPoint,1);
t21 = cell(nPoint,1);
t22 = cell(nPoint,1);
t11_R = cell(nPoint,1);
t12_R = cell(nPoint,1);
t21_R = cell(nPoint,1);
t22_R = cell(nPoint,1);
u11 = cell(nPoint,1);
u12 = cell(nPoint,1);
u21 = cell(nPoint,1);
u22 = cell(nPoint,1);

for loop1 = 1:nPoint
    
    % 非对角元素
    % loop1    
    r = [mat_r{loop1}];
    r_n = [mat_dr_D_dn{loop1}];
    r_x1 = [mat_r_x1{loop1}];
    r_x2 = [mat_r_x2{loop1}];
    
    u2=(-besselh(2,1,k2*r)+(k1/k2)^2*besselh(2,1,k1*r));
    u1_r=-(k2)*besselh(1,1,k2*r)+(2/k2)*(besselh(1,1,k2*r)-(k1/k2)*besselh(1,1,k1*r))./(r.*r)+((k1/k2)^2*besselh(0,1,k1*r)-besselh(0,1,k2*r))./r;
    u2_r=2*u1_r+besselh(1,1,k2*r)*k2+((k1^3)/(k2^2))*besselh(1,1,k1*r);

    t11_ji=1i/4*(((r_n+n1.*r_x1)+(a/b)*n1.*r_x1).*u1_r-((r_n+n1.*r_x1)+2*(n1.*r_x1-2*r_x1.*r_x1.*r_n)+(a/b)*n1.*r_x1).*(u2./r)-(2*r_x1.*r_x1.*r_n+(a/b)*n1.*r_x1).*u2_r);
    t12_ji=1i/4*((n1.*r_x2+(a/b)*n2.*r_x1).*u1_r-(n1.*r_x2+2*(n2.*r_x1-2*r_x1.*r_x2.*r_n)+(a/b)*n2.*r_x1).*(u2./r)-(2*r_x1.*r_x2.*r_n+(a/b)*n2.*r_x1).*u2_r);
    t21_ji=1i/4*((n2.*r_x1+(a/b)*n1.*r_x2).*u1_r-(n2.*r_x1+2*(n1.*r_x2-2*r_x1.*r_x2.*r_n)+(a/b)*n1.*r_x2).*(u2./r)-(2*r_x1.*r_x2.*r_n+(a/b)*n1.*r_x2).*u2_r);
    t22_ji=1i/4*(((r_n+n2.*r_x2)+(a/b)*n2.*r_x2).*u1_r-((r_n+n2.*r_x2)+2*(n2.*r_x2-2*r_x2.*r_x2.*r_n)+(a/b)*n2.*r_x2).*(u2./r)-(2*r_x2.*r_x2.*r_n+(a/b)*n2.*r_x2).*u2_r);

    t11_ji=arrW*t11_ji.*(0.5*ElemLen.');
    t11_ji(:,loop1)=0.5;
    t12_ji=arrW*t12_ji.*(0.5*ElemLen.');
    t12_ji(:,loop1)=0;
    t21_ji=arrW*t21_ji.*(0.5*ElemLen.');
    t21_ji(:,loop1)=0;
    t22_ji=arrW*t22_ji.*(0.5*ElemLen.');
    t22_ji(:,loop1)=0.5;

    t11_ji_R=t11_ji;
    t11_ji_R(loop1)=0;
    t12_ji_R=t12_ji;
    t12_ji_R(loop1)=0;
    t21_ji_R=t21_ji;
    t21_ji_R(loop1)=0;
    t22_ji_R=t22_ji;
    t22_ji_R(loop1)=0;
    
    t11(loop1) = {t11_ji};
    t12(loop1) = {t12_ji};
    t21(loop1) = {t21_ji};
    t22(loop1) = {t22_ji};

    t11_R(loop1) = {t11_ji_R};
    t12_R(loop1) = {t12_ji_R};
    t21_R(loop1) = {t21_ji_R};
    t22_R(loop1) = {t22_ji_R};

    u1=besselh(0,1,k2*r)-(1/k2)*besselh(1,1,k2*r)./r+(k1/k2)^2*(1/k1)*besselh(1,1,k1*r)./r;
    u2=-besselh(2,1,k2*r)+(k1/k2)^2*besselh(2,1,k1*r);
    u11_ji=1i/(4*b)*(u1-u2.*r_x1.*r_x1);
    u12_ji=1i/(4*b)*(-u2.*r_x1.*r_x2);
    u21_ji=1i/(4*b)*(-u2.*r_x1.*r_x2);
    u22_ji=1i/(4*b)*(u1-u2.*r_x2.*r_x2);

    u11_ji=arrW*u11_ji.*(0.5*ElemLen.');
    u11_ji(:,loop1)=0;
    u12_ji=arrW*u12_ji.*(0.5*ElemLen.');
    u12_ji(:,loop1)=0;
    u21_ji=arrW*u21_ji.*(0.5*ElemLen.');
    u21_ji(:,loop1)=0;
    u22_ji=arrW*u22_ji.*(0.5*ElemLen.');
    u22_ji(:,loop1)=0;

    % 对角元素  
   
    Lelem = ElemLen(loop1);
    r = Lelem*arrX.'/2;
    r_x1=-NormalVector(loop1,2);
    r_x2=NormalVector(loop1,1);
    A1=1+2*i*(log(r)+oula)/pi;
    B1=2*i/pi*(log(k2/2)-py_n);
    C1=2*i/pi*(log(k1/2)-py_n);
    D1=2*i/pi*(log(k2/2)-(k1/k2)^(2*n+2)*log(k1/2));
    E1=i/pi*(py_n+py_n1);
    u1_o=(i*(a+3*b)/(pi*(a+2*b)))*log(r);
    u2_o=i*(a+b)/(pi*(a+2*b));
    u1_d1=1/(2*(a+2*b))*(a+3*b+i/pi*(2*(a+2*b)*log(k2/2)+a+b+2*(a+3*b)*oula+2*b*log(k2/2)));
    u1_dd2=(-1)^x/(factorial(x))^2*(1/(2*(x+1))*(2*x+1+(k1/k2)^(2*x+2))*(k2*r/2).^(2*x).*A1+B1-1/(2*(x+1))*D1+1/(2*(x+1))*(1-(k1/k2)^(2*x+2))*E1);
    u2_dd1=(-1)^x/(factorial(x))^2*(x/(x+1)*(1-(k1/k2)^(2*x*2*x))*(k2*r/2).^(2*x).*A1+B1-(k1/k2)^(2*x+2)*C1-x/(x+1)*D1+1/(x+1)*(1-(k1/k2)^(2*x+2))*E1);
    u1_d2=0;
    u2_d=0;
for x=1:n
    u1_d2=u1_d2+u1_dd2;
    u2_d=u2_d+u2_dd1;
end
    u1_d=u1_d1+u1_d2;

    u11_j=i/(4*b)*(u1_d-(u2_o+u2_d)*r_x1^2)*Lelem/2;
    u12_j=-i/(4*b)*Lelem/2*(u2_o+u2_d)*r_x1*r_x2;
    u21_j=-i/(4*b)*Lelem/2*(u2_o+u2_d)*r_x1*r_x2;
    u22_j=i/(4*b)*Lelem/2*(u1_d-(u2_o+u2_d)*r_x2^2);
    
    u11_ji(:,loop1)=arrW*u11_j+Lelem*(a+3*b)/(4*b*pi*(a+2*b))*(log(2/Lelem)+1);
    u12_ji(:,loop1)=arrW*u12_j;
    u21_ji(:,loop1)=arrW*u21_j;
    u22_ji(:,loop1)=arrW*u22_j+Lelem*(a+3*b)/(4*b*pi*(a+2*b))*(log(2/Lelem)+1);

    u11(loop1) = {u11_ji};
    u12(loop1) = {u12_ji};
    u21(loop1) = {u21_ji};
    u22(loop1) = {u22_ji};

end
t11 = cell2mat(t11);
t12 = cell2mat(t12);
t21 = cell2mat(t21);
t22 = cell2mat(t22);
t11_R = cell2mat(t11_R);
t12_R = cell2mat(t12_R);
t21_R = cell2mat(t21_R);
t22_R = cell2mat(t22_R);
u11 = cell2mat(u11);
u12 = cell2mat(u12);
u21 = cell2mat(u21);
u22 = cell2mat(u22);
T1=zeros(nPoint,2*nPoint);
T2=zeros(nPoint,2*nPoint);
T3=zeros(2*nPoint,2*nPoint);
T1_R=zeros(nPoint,2*nPoint);
T2_R=zeros(nPoint,2*nPoint);
T3_R=zeros(2*nPoint,2*nPoint);
U1=zeros(nPoint,2*nPoint);
U2=zeros(nPoint,2*nPoint);
U3=zeros(2*nPoint,2*nPoint);
ii = 1:nPoint;
T1(:,2*ii-1)=t11(:,ii);
T1(:,2*ii)=t12(:,ii);
T2(:,2*ii-1)=t21(:,ii);
T2(:,2*ii)=t22(:,ii);
T3(2*ii-1,:)=T1(ii,:);
T3(2*ii,:)=T2(ii,:);

T1_R(:,2*ii-1)=t11_R(:,ii);
T1_R(:,2*ii)=t12_R(:,ii);
T2_R(:,2*ii-1)=t21_R(:,ii);
T2_R(:,2*ii)=t22_R(:,ii);
T3_R(2*ii-1,:)=T1_R(ii,:);
T3_R(2*ii,:)=T2_R(ii,:);

U1(:,2*ii-1)=u11(:,ii);
U1(:,2*ii)=u12(:,ii);
U2(:,2*ii-1)=u21(:,ii);
U2(:,2*ii)=u22(:,ii);
U3(2*ii-1,:)=U1(ii,:);
U3(2*ii,:)=U2(ii,:);

T=T3;
T_R=T3_R;
G=U3;