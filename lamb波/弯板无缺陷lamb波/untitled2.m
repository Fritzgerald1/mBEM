X = ff(1:2:length(ff))/1e5;

figure(1)

plot(X,abs(real(Rtot(1,:))),'r*-')
hold on
plot(X,-abs(imag(Rtot(1,:))),'r*--')
hold off
title('S0')
ylim([-0.03,0.03])

figure(2)

plot(X,real(Rtot(3,:)),'b*-')
hold on
plot(X,-imag(Rtot(3,:)),'b*--')
hold off
title('A0')

%{
nPropMode = 1:1:2*num;
figure('Name','散射系数')
plot(nPropMode,real(R_(nPropMode)),'r.-')
title('散射系数')
hold on
plot(nPropMode,imag(R_(nPropMode)),'r--')
hold off
xticks([1 3]); xticklabels({'S_0','A_0'})
legend('real','image')

nPropMode = 2*num+1:2:4*num;
figure('Name','透射系数')
plot(nPropMode,real(R_(nPropMode)),'r.-')
title('透射系数')
hold on
plot(nPropMode,imag(R_(nPropMode)),'r--')
hold off
xticks([1 3]); xticklabels({'S_0','A_0'})
legend('real','image')
%}