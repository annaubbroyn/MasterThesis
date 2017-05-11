option = 5;
lambda = [0.5 0.7 2];

system = 'SNS';
alpha = 0;%*pi/4;%0:pi/10:pi;
dalpha = 0;
alphaL = alpha;
alphaR = alpha + dalpha;

dy = 0.1;
dtheta = 0.1;
nB = 100;
W = 10;
x = 0;

ny = W/dy;
ntheta = lower(pi/dtheta);
n = nB;

if(option == 1)
    lambda = 1;
end

if(option == 1 || option == 5)
    B = linspace(-50/W,50/W,nB);
    %B = logspace(1,4,nB)/W;
else
    B = linspace(-50,50,nB);
end

figure(1)
hold on
for lam = 1:length(lambda)
    I = zeros(n,1);
    for k=1:n
        k
        y = -W/2;
        for i = 1:ny
            theta = -pi/2+dtheta;
            I_temp = zeros(n,1);
            for j = 1:ntheta-2                
                dI = localCurrent(x,y,theta,1/sqrt(B(k)),pi/2,option,lambda(lam),system,alphaL,alphaR);
                I(k,1) = I(k,1) + (2/W)*cos(theta)*dI*dy*dtheta;
                theta = theta + dtheta;
            end
            y = y + dy;
        end
    end
    if(option == 1 || option == 5)
        plot(B*W,abs(I(:,1)),'linewidth',1.5);    
        %plot(log2(B*W),log10(abs(I(:,1))),'linewidth',1.5);
    else
        plot(B,abs(I(:,1)),'linewidth',1.5);
    end
end

%Y-label
ylabel('$I_c/I_{c,0}$','interpret','latex','fontsize',20)
set(gca,'YTick',[0 0.5 1],'fontsize',20)

%X-label
if(option == 1 || option == 5)
    xlabel('$LW/l_m^2$','interpret','latex','fontsize',20)
    axis([B(1)*W B(length(B))*W 0 1]);
    set(gca,'XTick',[B(1)*W B(1)*W/2 0 -B(1)*W/2 -B(1)*W])
else
    xlabel('$L^2/l_m^2$','interpret','latex','fontsize',20)
    axis([B(1) B(length(B)) 0 1]);
    set(gca,'XTick',[B(1) B(1)/2 0 -B(1)/2 -B(1)])
end

%Legend
if(option ~= 1)
    legendInfo = cell(length(lambda),1);
    for i=1:length(lambda)
        lambdaString = [num2str(floor(lambda(i))) '.' num2str(floor(10*(lambda(i)-floor(lambda(i)))))];
        legendInfo{i} = ['$\lambda/L = $' lambdaString];
    end
    h = legend(legendInfo);
    set(h,'Interpreter','Latex','fontsize',12,'position',[0.69 0.725 0.2 0.2])
end

filename = ['C:\Users\Anna\Documents\GitHub\MasterThesis\Matlab\Figures\' system '\Dist' num2str(option) '\CriticalCurrent\Critical' num2str(option) 'nB' num2str(nB) 'dy' strrep(num2str(dy),'.','-') 'dtheta' strrep(num2str(dtheta),'.','-')];
print(filename,'-dpng')
savefig(filename)

%For log-plot:
%set(gca,'YTick',[-3 -2 -1 0],'YtickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^0'},'fontsize',20)
%set(gca,'XTick',[log2(10) log2(20) log2(30) log2(40)],'XtickLabel',{'10','20','30','40'},'fontsize',20)
%axis([log2(10) log2(40) -2 0]);
