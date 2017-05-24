l= 0.5;
option = 1;
lambda = 2;%0.1:0.1:1;
video = 0;
vortex = 0;
system = 'Dwave';
alpha = 0;%0:pi/10:pi;
dalpha = 0;
alphaL = alpha;
alphaR = alpha+dalpha;
swave = 0;
zeeman = 0.0075;

phi = 0;
nx = 50;
ny = 400;
%nx = 5;
%ny = 40;

ntheta = 100;
dimy = 2;
l_arrow = 0.2;

cm = colormap(jet);
cm = cm(10:50,:);

if(video == 1)
    writerobj = VideoWriter('out.avi');
    writerobj.FrameRate = 5;
    open(writerobj);
end

if(option == 1)
    lambda = 1;
end
if( ~strcmp(system,'Dwave'))
    alphaL = 0;
    alphaR = 0;
end

for k_alpha = 1:length(alpha)
    %toPrint=['alpha: ',num2str(k_alpha),'/',num2str(length(alpha))];
    %disp(toPrint)

    for k = 1:length(lambda)

        %toPrint=['lambda: ',num2str(k),'/',num2str(length(lambda))];
        %disp(toPrint)

        dx = 1/nx;
        dy = 2*dimy/ny;
        dv = 0.25*pi*l^2;
        if(dv>4)
            dv = 0.1;
        end
        nv0 = ceil(2/pi*(phi+2*dimy/l^2));
        yv0=0.5*l^2*(phi-nv0*pi/2);
        while(dv>0 && yv0<-dimy)
            yv0 = yv0 + 0.25*dv;
        end
        xv0 = -dv*floor(0.5/dv);
        ny0 = floor((yv0+dimy)/dy)+1;
        nx0 = floor((xv0+0.5)/dx)+1;
        dnx = floor(dv/dx)+1;
        dny = floor(dv/dy)+1;
        dnx1 = 1.1*floor(dnx);
        dny1 = floor(dny);
        nx1 = nx;
        ny1 = ny;

        x0 = linspace(-0.5,0.5,nx);
        y0 = linspace(-dimy,dimy,ny);
        jx = zeros(nx,ny);
        jy = zeros(nx,ny);
        jtot = zeros(nx,ny);
        theta = linspace(-pi/2,pi/2,ntheta);

        [y,x] = meshgrid(y0,x0);

        I = linspace(0,0,nx);

        for i=1:nx
            %disp([num2str(i) '/' num2str(nx)])
            for j = 1:ny
                %disp([num2str(j) '/' num2str(ny)])
                fun= @(theta)localCurrent(x(i,j),y(i,j),theta, l, phi,option,lambda(k),system,alphaL(k_alpha),alphaR(k_alpha),zeeman,swave);
                dtheta = pi/ntheta;
                theta = -pi/2 + dtheta;
                for(n = 1:98)
                    jx(i,j) = jx(i,j)+fun(theta)*cos(theta)*dtheta;
                    jy(i,j) = jy(i,j)+fun(theta)*sin(theta)*dtheta;
                    theta = theta + dtheta;
                end

                jtot(i,j) = sqrt(jx(i,j)^2+jy(i,j)^2);
            end
        end

        fig = figure(k);
        hold on
        fig.Units = 'normalized';
        axis([-0.5 0.5 -dimy dimy+.1]);
        fig.OuterPosition = [.05 .05 .1 .9];
        L = sqrt(jx.^2 + jy.^2);
        
%         quiver(x,y,jx,jy,'r','AutoScaleFactor',.3);
%         quiver(x(1:2:nx,1:2:ny),y(1:2:nx,1:2:ny),jx(1:2:nx,1:2:ny),jy(1:2:nx,1:2:ny),'r','AutoScaleFactor',0.3);
%         quiver(x(1:3:nx,1:3:ny),y(1:3:nx,1:3:ny),jx(1:3:nx,1:3:ny),jy(1:3:nx,1:3:ny),'r','AutoScaleFactor',0.4);
         quiver(x(1:3:nx,1:3:ny),y(1:3:nx,1:3:ny),jx(1:3:nx,1:3:ny),jy(1:3:nx,1:3:ny),'r','AutoScaleFactor',0.3);
%         quiver(x(1:5:nx,1:5:ny),y(1:5:nx,1:5:ny),jx(1:5:nx,1:5:ny),jy(1:5:nx,1:5:ny),'r','AutoScaleFactor',0.5);
%         quiver(x(1:4:nx,1:4:ny),y(1:4:nx,1:4:ny),jx(1:4:nx,1:4:ny)./L(1:4:nx,1:4:ny),jy(1:4:nx,1:4:ny)./L(1:4:nx,1:4:ny),'r','AutoScaleFactor',.4);
%         quiver(x(1:4:nx,1:10:ny),y(1:4:nx,1:10:ny),jx(1:4:nx,1:10:ny)./L(1:4:nx,1:10:ny),jy(1:4:nx,1:10:ny)./L(1:4:nx,1:10:ny),'r','AutoScaleFactor',.5);
%          jx(abs(jx)>0.3) = 0;
%          jy(abs(jy)>0.3) = 0;
%              if(dv>0)
%                  quiver(x(nx0:dnx1:nx1,ny0:dny1:ny1),y(nx0:dnx1:nx1,ny0:dny1:ny1),jx(nx0:dnx1:nx1,ny0:dny1:ny1),jy(nx0:dnx1:nx1,ny0:dny1:ny1),'r','AutoScaleFactor',l_arrow);
%              end
%         
        h = surf(x,y,0.5*jtot);
        view(0,90); shading interp
        set(h,'ZData',jtot-1);
        colormap(cm)
        hcb = colorbar('southoutside');
        set(hcb,'YTick',[-1 -0.75 -0.5],'YtickLabel',{'0','0.25','0.5'},'fontsize',10)
        set(gca,'XTick',[-0.5 0 0.5],'XtickLabel',{'-L/2','0','L/2'},'fontsize',10);
        set(gca,'YTick',[-2 -1 0 1 2],'YtickLabel',{'-2L','-L','0','L','2L'},'fontsize',10);
        caxis([-1 -0.5])
        handle = title(['$\lambda/L = $' num2str(lambda(k))],'interpret','latex');
        set(handle,'Units','Normalized','Position',[0.5 0.98 1]);

        figGcf = gcf;
        figGcf.PaperUnits = 'normalized';
        figGcf.PaperPosition = [0 0 .4 1];

        if(video == 0)
            filename1 = ['C:\Users\Anna\Documents\GitHub\MasterThesis\Matlab\Figures\' system '\Dist' num2str(option) '\Dist' num2str(option) '_l_0-' num2str(l*10) '_lambda_' num2str(floor(lambda(k))) '-' num2str(floor(100*(lambda(k)-floor(lambda(k))))) 'phi_pi-' floor(num2str(pi/phi)) '_alphaL-' strrep(num2str(alphaL),'.','-') '_alphaR-' strrep(num2str(alphaR),'.','-') '_h_' strrep(num2str(zeeman),'.','-') '_swave_' strrep(num2str(swave),'.','-')]; 
            %filename2 = ['C:\Users\Anna\Documents\GitHub\MasterThesis\Matlab\Figures\Dist' num2str(option) '\Matlab\Dist' num2str(option) '_l_0-' num2str(l*10) '_lambda_' num2str(floor(lambda(k))) '-' num2str(floor(100*(lambda(k)-floor(lambda(k))))) 'phi_pi-' floor(num2str(pi/phi))]; 
            disp(filename1)
            print(fig,filename1,'-dpng');
            %print(filename1,'-dpng')
            %savefig(filename2)
        else
            frame = getframe(figGcf);
            writeVideo(writerobj,frame);
        end

        close(fig);
        if(vortex == 1)
            hold off

            fig2 = figure(2*k);
            hold on
            jx(jx<0.05) = 0;
            jx(jx>0.3) = 0;
            jy(jy<0.05) = 0;
            jy(jy>0.3) = 0;

            jx(abs(jy)>abs(jx)) = 0;
            jy(abs(jx)>abs(jy)) = 0;
            jx(abs(jx)<0.1) = 0;
            jy(abs(jy)<0.05) = 0;
            jx(abs(jx)>0.2) = 0;
            jy(abs(jy)>0.2) = 0;
            L = sqrt(jx.^2 + jy.^2);
% 
            quiver(x(1:6:nx,1:10:ny),y(1:6:nx,1:10:ny),jx(1:6:nx,1:10:ny)./L(1:6:nx,1:10:ny),jy(1:6:nx,1:10:ny)./L(1:6:nx,1:10:ny),'k','AutoScaleFactor',.5,'AlignVertexCenters','on','linewidth',1);
%             quiver(x(1:8:nx,1:8:ny),y(1:8:nx,1:8:ny),jx(1:8:nx,1:8:ny)./L(1:8:nx,1:8:ny),jy(1:8:nx,1:8:ny)./L(1:8:nx,1:8:ny),'k','AutoScaleFactor',.4,'linewidth',1);
%             quiver(x(1:8:nx,1:8:ny),y(1:8:nx,1:8:ny),jx(1:8:nx,1:8:ny),jy(1:8:nx,1:8:ny),'k','AutoScaleFactor',.4,'linewidth',1);
%             quiver(x(1:3:nx,1:3:ny),y(1:3:nx,1:3:ny),jx(1:3:nx,1:3:ny),jy(1:3:nx,1:3:ny),'k','AutoScaleFactor',.3,'linewidth',1);
            

            dnx2 = floor(1.25*dnx);
            dny2 = floor(1.25*dny);
            
            
%             quiver(x(nx0:dnx2:nx1,ny0:dny2:ny1),y(nx0:dnx2:nx1,ny0:dny2:ny1),jx(nx0:dnx2:nx1,ny0:dny2:ny1),jy(nx0:dnx2:nx1,ny0:dny2:ny1),'k','AutoScaleFactor',.3,'linewidth',1);
%             quiver(x(1:3:nx,1:2:ny),y(1:3:nx,1:2:ny),jx(1:3:nx,1:2:ny)./L(1:3:nx,1:2:ny),jy(1:3:nx,1:2:ny)./L(1:3:nx,1:2:ny),'k','AutoScaleFactor',.4,'linewidth',1);
% 
%             plot(0,-pi*l^2*1.25,'ob','markerfacecolor','b','markersize',10)
%             plot(0,-pi*l^2*0.75,'or','markerfacecolor','r','markersize',10)
%             plot(0,-pi*l^2*0.25,'ob','markerfacecolor','b','markersize',10)
%             plot(0,pi*l^2*0.25,'or','markerfacecolor','r','markersize',10)
%             plot(0,pi*l^2*0.75,'ob','markerfacecolor','b','markersize',10)
%             plot(0,pi*l^2*1.25,'or','markerfacecolor','r','markersize',10)
    

            fig2.Units = 'normalized';
            fig2.GraphicsSmoothing = 'on';
            fig2.Renderer = 'opengl';
            axis([-0.5 0.5 -dimy dimy+.1]);
            fig2.OuterPosition = [.05 .05 .17 .9];
            set(gca,'XTick',[-0.5 0 0.5],'XtickLabel',{'-L/2','0','L/2'},'fontsize',15);
            set(gca,'YTick',[-2 -1 0 1 2],'YtickLabel',{'-2L','-L','0','L','2L'},'fontsize',15);
            handle = title(['$\lambda/L = $' num2str(lambda(k))],'interpret','latex');
            set(handle,'Units','Normalized','Position',[0.5 0.98 1]);

            figGcf = gcf;
            figGcf.PaperUnits = 'normalized';
            figGcf.PaperPosition = [0 0 .4 1];

            filename1 = ['C:\Users\Anna\Documents\GitHub\MasterThesis\Matlab\Figures\' system '\Dist' num2str(option) '\Vortex' num2str(option) '_l_0-' num2str(l*10) '_lambda_' num2str(floor(lambda(k))) '-' num2str(floor(100*(lambda(k)-floor(lambda(k))))) 'phi_pi-' floor(num2str(pi/phi)) '_alphaL-' strrep(num2str(alphaL),'.','-') '_alphaR-' strrep(num2str(alphaR),'.','-') '_h_' strrep(num2str(zeeman),'.','-') '_swave_' strrep(num2str(swave),'.','-')]; 
            %filename1 = ['C:\Users\Anna\Documents\GitHub\Project\Figures\Dist' num2str(option) '\Vortex' num2str(option) '_l_0-' num2str(l*10) '_lambda_' num2str(floor(lambda(k))) '-' num2str(floor(100*(lambda(k)-floor(lambda(k))))) 'phi_pi-' floor(num2str(pi/phi))]; 
            %filename2 = ['C:\Users\Anna\Documents\GitHub\Project\Figures\Dist' num2str(option) '\Matlab\Vortex' num2str(option) '_l_0-' num2str(l*10) '_lambda_' num2str(floor(lambda(k))) '-' num2str(floor(100*(lambda(k)-floor(lambda(k))))) 'phi_pi-' floor(num2str(pi/phi))]; 
            print(filename1,'-dpng')
            %print(filename1,'-dpdf')
            %savefig(filename2)
        end

    end
end
if(video == 1)
    close(writerobj)
end