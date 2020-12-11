clear all; close all; clc;
load wind;

Nx = 128;
Ny = 64;

x = linspace(0,2,Nx+1);
y = linspace(0,1,Ny+1);
xm = x(1:end-1)+(x(2)-x(1))/2;
ym = y(1:end-1)+(y(2)-y(1))/2;

length = 2051;

for i = 2:length
    pdata = readmatrix(append('post/niubility',num2str(i,'%d'),'.dat'));
    edata = readmatrix(append('post/ensemble',num2str(i,'%d'),'.dat'));
    
    x = reshape(pdata(:,1),[Nx,Ny]);
    y = reshape(pdata(:,2),[Nx,Ny]);

    u = reshape(pdata(:,3),[Nx,Ny]);
    v = reshape(pdata(:,4),[Nx,Ny]);
    p = reshape(pdata(:,5),[Nx,Ny]);
    
    xp = edata(:,1);
    yp = edata(:,2);
    up = edata(:,3);
    vp = edata(:,4);
    
    figure(1)
    pcolor(x,y,sqrt(u.^2+v.^2))
    shading interp
    hold on
    scatter(xp,yp,sqrt(up.^2+vp.^2),'MarkerEdgeColor','k')
%     contourf(x,y,u)
%     hs = streamslice(x',y',u',v');
%     set(hs,'color','w','linewidth',1);
    set(gcf,'Color',[1,1,1]);
    set(gcf,'Position',[200 500 2000 400])
    c = colorbar;
    set(c,'Limits',[0,0.1])
    ylim([0.3 0.7])
    
%     caxis([0 0.0001])
    hold off
end
