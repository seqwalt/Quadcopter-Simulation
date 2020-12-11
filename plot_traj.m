function plot_traj(x,z,th,snaps,xlab,ylab,titl)

hold on
axis equal
plot(x,z);
x_range = abs(max(x)-min(x));
z_range = abs(max(z)-min(z));
traj_len = sqrt(x_range^2 + z_range^2);
even_time = 1:int16(size(th,1)/snaps):size(th,1);

for i = even_time
    R = [cos(-th(i)), -sin(-th(i));
         sin(-th(i)), cos(-th(i))];
    L = R*[-.03*traj_len;0] + [x(i);z(i)];
    R = R*[ .03*traj_len;0] + [x(i);z(i)];
    quad = [L,R];
    if i == even_time(1)
       plot(quad(1,:),quad(2,:),'g','LineWidth',4);
    elseif i == even_time(end)
        plot(quad(1,:),quad(2,:),'r','LineWidth',4);
    else
        plot(quad(1,:),quad(2,:),'b','LineWidth',3);
    end
end
ax = gca;
ax.FontSize = 14;
xl = xlabel(xlab,'Interpreter','latex'); xl.FontSize = 18;
yl = ylabel(ylab,'Interpreter','latex'); yl.FontSize = 18;
ti = title(titl,'Interpreter','latex');  ti.FontSize = 18;