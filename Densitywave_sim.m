
clear
close all
n = 15000;     
nt = 2000;     

v=5;          
y=zeros(n,nt);    
x=zeros(n,nt);    
vy = zeros(n,nt);
vx = zeros(n,nt);
theta = zeros(n,nt);

%% 
load('att_t.mat')
NPD=lognfit(Y);
mu = NPD(1);
sigma = NPD(2);
att1 = lognrnd(mu,sigma,1,n);
att=round(att1);

lfx = [1:1:60];
lfy = zeros(1,60);
for i = 1:60
    lfy(i) = 1/(lfx(i)*sigma*sqrt(2*pi))*exp(-(log(lfx(i))-mu)^2/(2*sigma*sigma));
end
wt=15*rand(1,n)+5*att;
wt = round(wt);

%% 
Nc = 100;
Pn = zeros(n,Nc);
Nn = zeros(n,Nc);

load('swim_t.mat')

SimPara(4,i) = 0.15 * SimPara(2,i);
SimPara(5,i) = 0.15 * SimPara(3,i);

for i = 1:Nc
    Pn(:,i) = normrnd(SimPara(2,i),SimPara(4,i),n,1);
    Pn(:,i) = round(abs(Pn(:,i)));

    Nn(:,i) = normrnd(SimPara(3,i),SimPara(5,i),n,1);
    Nn(:,i) = round(abs(Nn(:,i)));
end

%%
Hb = 150;        
Lb = -150;
Tb = -60;
Rb = 60;

y(:,1)=(rand(n,1)-0.5)*2*Hb;    
x(:,1)=(rand(n,1)-0.5)*2*Rb;
theta(:,1)=rand(n,1)*2*pi;

vx(:,1)=v*cos(theta(:,1));
vy(:,1)=v*sin(theta(:,1));


%% 
kb = 1.380649*10^(-23);
T = 298.15;
mu = 8.9*10^(-4);
R = 10*10^(-6);
Dt = kb*T/(6*pi*mu*R)*10^(6)*10^(6);
Dr = kb*T/(8*pi*mu*R*R*R);

dx = sqrt(2*Dt)*wgn(n,nt,1);
dy = sqrt(2*Dt)*wgn(n,nt,1);
dr = sqrt(2*Dr)*wgn(n,nt,1);

dtheta1 = wgn(n,200,1)/2;
dtheta2 = wgn(n,nt-200,1)/2;
dtheta = [dtheta1,dtheta2];

%%
Ap1 = rand(n,Nc);
Ap2 = rand(n,Nc);
Cp =1;

%%
for i = 1:n
    A = 2;
    B = Pn(i,1) + 1;
    for w = 1:(Nc-1)
        z=0;
        for u = A:B   
            if u>nt
                break;
            end
%% 
            if Ap1(i,w)<=Cp && u==(A+att(i))
                for q = u:(u+wt(i))
                    x(i,q) = x(i,q-1);
                    y(i,q) = y(i,q-1);
                end
                for p = (u+wt(i)+1):(wt(i)+1+B)

                    if p>nt
                        break;
                    end

                    k = 0.5;
                    thetaI = 3*pi/2;
                    vx(i,p)=v*cos(theta(i,p-1)) + dx(i,p-1);     
                    vy(i,p)=v*sin(theta(i,p-1)) + dy(i,p-1);
                    theta(i,p)= theta(i,p-1) + k*(thetaI-theta(i,p-1)) + dtheta(i,p-1) + dr(i,p-1);

                    x(i,p) = x(i,p-1) + vx(i,p);
                    y(i,p) = y(i,p-1) + vy(i,p);

                    if x(i,p)>Rb
                        x(i,p)=Rb;
                    end
                    if x(i,p)<Tb
                        x(i,p)=Tb;
                    end

                    if y(i,p)<Lb
                        y(i,p)=Lb;
                    end
                    if y(i,p)>Hb
                        y(i,p)=Hb;
                    end
                end
                z=wt(i);
                break
            end

            k = 0.5;
            thetaI = 3*pi/2;
            vx(i,u)=v*cos(theta(i,u-1)) + dx(i,u-1);     
            vy(i,u)=v*sin(theta(i,u-1)) + dy(i,u-1);
            theta(i,u)= theta(i,u-1) + k*(thetaI-theta(i,u-1)) + dtheta(i,u-1) + dr(i,u-1);

            x(i,u) = x(i,u-1) + vx(i,u);
            y(i,u) = y(i,u-1) + vy(i,u);

            if x(i,u)>Rb
                 x(i,u)=Rb;
            end
            if x(i,u)<Tb
                 x(i,u)=Tb;
            end
            if y(i,u)<Lb
                y(i,u)=Lb;
            end
            if y(i,u)>Hb
                y(i,u)=Hb;
            end

        end
            A = B+1+z;
            B = B+Nn(i,w)+z;
            z=0;
        
        for u = A:B
            if u>nt
                break;
            end

            if u>nt
                break;
            end

            k = 0.5;
            thetaI = 1*pi/2;
            vx(i,u)=v*cos(theta(i,u-1)) + dx(i,u-1);      
            vy(i,u)=v*sin(theta(i,u-1)) + dy(i,u-1);
            theta(i,u)= theta(i,u-1) + k*(thetaI-theta(i,u-1)) + dtheta(i,u-1) + dr(i,u-1);

            x(i,u) = x(i,u-1) + vx(i,u);
            y(i,u) = y(i,u-1) + vy(i,u);

            if x(i,u)>Rb
                 x(i,u)=Rb;
            end
            if x(i,u)<Tb
                 x(i,u)=Tb;
            end
            if y(i,u)<Lb
                y(i,u)=Lb;
            end
            if y(i,u)>Hb
                y(i,u)=Hb;
            end

        end
            A = B+1+z;
            B = B+Pn(i,w+1)+z;
            z=0;
        
    end
end


%% 
AA4 = zeros(61,60);
v = 0;
ddy = 5;
dn = (Hb-Lb)/ddy;
nn = n;
Ln = y;
dmap = zeros(dn,nt); 

for u = 1:30:1801
    for i = 1:dn
          if  i==dn
              k = find(  Ln(:,u) >= Lb+(i-1)*ddy &  Ln(:,u) <= Lb+i*ddy );
          else
              k = find(  Ln(:,u) >= Lb+(i-1)*ddy &  Ln(:,u) < Lb+i*ddy );
          end
        dmap(i,u) = size(k,1)/nn ;
    end
end

for u = 1:30:1801
        v = v+1;
        AA4(v,:) = dmap(:,u);
end

figure(1)
colormap(jet(256))
h1 = imagesc(AA4,'AlphaData',~isnan(AA4));
caxis([0, 0.03]);
c1 = colorbar;
set(c1,'XTick',0:0.01:0.03,'fontsize',21);
set(c1,'XTickLabel',{'0','0.01','0.02','0.03'})
xlabel('\it{x (mm)}');
ylabel('\it{t (mins)}');
xticks(1:10:61)
yticks(1:10:61)
set(gca,'XTickLabel',{'0','','1','','2','','3'})
set(gca,'YTickLabel',{'0','','10','','20','','30'})
set(gca,'TickDir','out');
set(gca,'linewidth',1.2)
set(gca,'FontWeight','normal')
set(gca,'FontSize',21)
set(gca,'XColor',[0,0,0])
set(gca,'YColor',[0,0,0])
set(gca,'ZColor',[0,0,0])
hold on
