% (BA)^{N}(BCB) ^{1}(AB)^{N}ÐÍ
clear all;k1=2,k2=3,k3=4,k4=5,;n=1;
n0=1;na=3.230;nb=1.35;nc=1.0;
c1=0/6; %½Ç¶È¿Éµ÷
ca=asin(n0*sin(c1)/na);cb=asin(na*sin(ca)/nb);cc=asin(nb*sin(cb)/nc);
c6=c1;

d0=87.5;da=1480/(4*na);db=1480/(4*nb);dc=750;d4=87.5; %µ¥Î»nm
p1=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n0*cos(c1); %¿ÕÆø²ã
pa=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*na*cos(ca); %A²ã
pb=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nb*cos(cb); %B²ã
pc=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nc*cos(cc); %C²ã
p4=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n0*cos(c6); %ÓÒ¿ÕÆø²ã

t1=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc=2*pi*nc*dc*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc),-i*sin(Fc)/pc;-i*pc*sin(Fc),cos(Fc)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

M= M1*((Mb*Ma)^k1*(Mb*Mc*Mb)^n*(Ma*Mb)^k1)*M4;
R=(abs(((M(1,1)+M(1,2)*p4)*p1-(M(2,1)+M(2,2)*p4))/((M(1,1)+M(1,2)*p4)*p1+(M(2,1)+M(2,2)*p4))))^2;
T=1-R;
t1=[t1,T];
end

t2=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc=2*pi*nc*dc*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc),-i*sin(Fc)/pc;-i*pc*sin(Fc),cos(Fc)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

M= M1*((Mb*Ma)^k2*(Mb*Mc*Mb)^n*(Ma*Mb)^k2)*M4;
R=(abs(((M(1,1)+M(1,2)*p4)*p1-(M(2,1)+M(2,2)*p4))/((M(1,1)+M(1,2)*p4)*p1+(M(2,1)+M(2,2)*p4))))^2;
T=1-R;
t2=[t2,T];
end

t3=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc=2*pi*nc*dc*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc),-i*sin(Fc)/pc;-i*pc*sin(Fc),cos(Fc)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

M= M1*((Mb*Ma)^k3*(Mb*Mc*Mb)^n*(Ma*Mb)^k3)*M4;
R=(abs(((M(1,1)+M(1,2)*p4)*p1-(M(2,1)+M(2,2)*p4))/((M(1,1)+M(1,2)*p4)*p1+(M(2,1)+M(2,2)*p4))))^2;
T=1-R;
t3=[t3,T];
end

t4=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc=2*pi*nc*dc*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc),-i*sin(Fc)/pc;-i*pc*sin(Fc),cos(Fc)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

M= M1*((Mb*Ma)^k4*(Mb*Mc*Mb)^n*(Ma*Mb)^k4)*M4;
R=(abs(((M(1,1)+M(1,2)*p4)*p1-(M(2,1)+M(2,2)*p4))/((M(1,1)+M(1,2)*p4)*p1+(M(2,1)+M(2,2)*p4))))^2;
T=1-R;
t4=[t4,T];
end

figure;
subplot(2,2,1)
plot(g,t1,'k')
title('N=2');
ylabel(' Transmission');
xlabel('¦Ë[nm]');
subplot(2,2,2)
plot(g,t2,'k')
title('N=3');
ylabel(' Transmission');
xlabel('¦Ë[nm]');
subplot(2,2,3)
plot(g,t3,'k')
title('N=4');
ylabel(' Transmission');
xlabel('¦Ë[nm]');
subplot(2,2,4)
plot(g,t4,'k')
title('N=5');
ylabel(' Transmission');
xlabel('¦Ë[nm]');

