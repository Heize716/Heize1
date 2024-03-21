% (BA)^{4}(BCB) ^{1}(AB)^{4}ÐÍ
clear all;k=4;n=1;
n0=1;na=3.230;nb=1.35;nc=1.0;
c1=0/6; %½Ç¶È¿Éµ÷
ca=asin(n0*sin(c1)/na);cb=asin(na*sin(ca)/nb);cc=asin(nb*sin(cb)/nc);
c6=c1;

d0=87.5;da=1480/(4*na);db=1480/(4*nb);dc1=600;dc2=650;dc3=700;dc4=750;dc5=800;dc6=850;d4=87.5; %µ¥Î»nm
%TM Mode
%p1=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n0/cos(c1); %¿ÕÆø²ã
%pa=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*na/cos(ca); %A²ã
%pb=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nb/cos(cb); %B²ã
%pc=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nc/cos(cc); %C²ã
%p4=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n/cos(c6); %ÓÒ¿ÕÆø²ã
%TE Mode
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
Fc=2*pi*nc*dc1*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc),-i*sin(Fc)/pc;-i*pc*sin(Fc),cos(Fc)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

M= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
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
Fc2=2*pi*nc*dc2*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc2),-i*sin(Fc2)/pc;-i*pc*sin(Fc2),cos(Fc2)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MB= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R2=(abs(((MB(1,1)+MB(1,2)*p4)*p1-(MB(2,1)+MB(2,2)*p4))/((MB(1,1)+MB(1,2)*p4)*p1+(MB(2,1)+MB(2,2)*p4))))^2;
T2=1-R2;
t2=[t2,T2];
end

t3=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc3=2*pi*nc*dc3*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc3),-i*sin(Fc3)/pc;-i*pc*sin(Fc3),cos(Fc3)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MC= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R3=(abs(((MC(1,1)+MC(1,2)*p4)*p1-(MC(2,1)+MC(2,2)*p4))/((MC(1,1)+MC(1,2)*p4)*p1+(MC(2,1)+MC(2,2)*p4))))^2;
T3=1-R3;
t3=[t3,T3];
end

t4=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc4=2*pi*nc*dc4*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc4),-i*sin(Fc4)/pc;-i*pc*sin(Fc4),cos(Fc4)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MD= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R4=(abs(((MD(1,1)+MD(1,2)*p4)*p1-(MD(2,1)+MD(2,2)*p4))/((MD(1,1)+MD(1,2)*p4)*p1+(MD(2,1)+MD(2,2)*p4))))^2;
T4=1-R4;
t4=[t4,T4];
end

t5=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc5=2*pi*nc*dc5*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc5),-i*sin(Fc5)/pc;-i*pc*sin(Fc5),cos(Fc5)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

ME= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R5=(abs(((ME(1,1)+ME(1,2)*p4)*p1-(ME(2,1)+ME(2,2)*p4))/((ME(1,1)+ME(1,2)*p4)*p1+(ME(2,1)+ME(2,2)*p4))))^2;
T5=1-R5;
t5=[t5,T5];
end

t6=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc6=2*pi*nc*dc6*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc6),-i*sin(Fc6)/pc;-i*pc*sin(Fc6),cos(Fc6)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MF= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R6=(abs(((MF(1,1)+MF(1,2)*p4)*p1-(MF(2,1)+MF(2,2)*p4))/((MF(1,1)+MF(1,2)*p4)*p1+(MF(2,1)+MF(2,2)*p4))))^2;
T6=1-R6;
t6=[t6,T6];
end

% (BA)^{5}(BCB) ^{1}(AB)^{5}ÐÍ
k=5;n=1;
n0=1;na=3.230;nb=1.35;nc=1.0;
c1=0/6; %½Ç¶È¿Éµ÷
ca=asin(n0*sin(c1)/na);cb=asin(na*sin(ca)/nb);cc=asin(nb*sin(cb)/nc);
c6=c1;

d0=87.5;da=1480/(4*na);db=1480/(4*nb);dc1=600;dc2=650;dc3=700;dc4=750;dc5=800;dc6=850;d4=87.5; %µ¥Î»nm
%TM Mode
%p1=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n0/cos(c1); %¿ÕÆø²ã
%pa=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*na/cos(ca); %A²ã
%pb=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nb/cos(cb); %B²ã
%pc=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nc/cos(cc); %C²ã
%p4=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n/cos(c6); %ÓÒ¿ÕÆø²ã
%TE Mode
p1=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n0*cos(c1); %¿ÕÆø²ã
pa=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*na*cos(ca); %A²ã
pb=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nb*cos(cb); %B²ã
pc=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*nc*cos(cc); %C²ã
p4=(8.854*10^(-12)/(4*pi*10^(-7)))^(1/2)*n0*cos(c6); %ÓÒ¿ÕÆø²ã

t_1=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc=2*pi*nc*dc1*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc),-i*sin(Fc)/pc;-i*pc*sin(Fc),cos(Fc)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

M= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R=(abs(((M(1,1)+M(1,2)*p4)*p1-(M(2,1)+M(2,2)*p4))/((M(1,1)+M(1,2)*p4)*p1+(M(2,1)+M(2,2)*p4))))^2;
T=1-R;
t_1=[t_1,T];
end

t_2=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc2=2*pi*nc*dc2*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc2),-i*sin(Fc2)/pc;-i*pc*sin(Fc2),cos(Fc2)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MB= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R2=(abs(((MB(1,1)+MB(1,2)*p4)*p1-(MB(2,1)+MB(2,2)*p4))/((MB(1,1)+MB(1,2)*p4)*p1+(MB(2,1)+MB(2,2)*p4))))^2;
T2=1-R2;
t_2=[t_2,T2];
end

t_3=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc3=2*pi*nc*dc3*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc3),-i*sin(Fc3)/pc;-i*pc*sin(Fc3),cos(Fc3)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MC= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R3=(abs(((MC(1,1)+MC(1,2)*p4)*p1-(MC(2,1)+MC(2,2)*p4))/((MC(1,1)+MC(1,2)*p4)*p1+(MC(2,1)+MC(2,2)*p4))))^2;
T3=1-R3;
t_3=[t_3,T3];
end

t_4=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc4=2*pi*nc*dc4*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc4),-i*sin(Fc4)/pc;-i*pc*sin(Fc4),cos(Fc4)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MD= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R4=(abs(((MD(1,1)+MD(1,2)*p4)*p1-(MD(2,1)+MD(2,2)*p4))/((MD(1,1)+MD(1,2)*p4)*p1+(MD(2,1)+MD(2,2)*p4))))^2;
T4=1-R4;
t_4=[t_4,T4];
end

t_5=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc5=2*pi*nc*dc5*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc5),-i*sin(Fc5)/pc;-i*pc*sin(Fc5),cos(Fc5)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

ME= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R5=(abs(((ME(1,1)+ME(1,2)*p4)*p1-(ME(2,1)+ME(2,2)*p4))/((ME(1,1)+ME(1,2)*p4)*p1+(ME(2,1)+ME(2,2)*p4))))^2;
T5=1-R5;
t_5=[t_5,T5];
end

t_6=[];
g=linspace(1000,2500,10000);%nm
for x=g
F1=2*pi*n0*d0*cos(c1)/x;
Fa=2*pi*na*da*cos(ca)/x;
Fb=2*pi*nb*db*cos(cb)/x;
Fc6=2*pi*nc*dc6*cos(cc)/x;
F4=2*pi*n0*d4*cos(c6)/x; 
M1=[cos(F1),-i*sin(F1)/p1;-i*p1*sin(F1),cos(F1)]; %×ó¿ÕÆø²ã
Ma=[cos(Fa),-i*sin(Fa)/pa;-i*pa*sin(Fa),cos(Fa)]; %A²ã
Mb=[cos(Fb),-i*sin(Fb)/pb;-i*pb*sin(Fb),cos(Fb)]; %B²ã
Mc=[cos(Fc6),-i*sin(Fc6)/pc;-i*pc*sin(Fc6),cos(Fc6)]; %C²ã
M4=[cos(F4),-i*sin(F4)/p4;-i*p4*sin(F4),cos(F4)]; %ÓÒ¿ÕÆø²ã

MF= M1*((Mb*Ma)^k*(Mb*Mc*Mb)^n*(Ma*Mb)^k)*M4;
R6=(abs(((MF(1,1)+MF(1,2)*p4)*p1-(MF(2,1)+MF(2,2)*p4))/((MF(1,1)+MF(1,2)*p4)*p1+(MF(2,1)+MF(2,2)*p4))))^2;
T6=1-R6;
t_6=[t_6,T6];
end

%Í¸ÉäÂÊ-²¨³¤
figure
subplot(1,2,1)
plot(g,t1,'k',g,t2,'r',g,t3,'g',g,t4,'b',g,t5,'y',g,t6,'m')
legend('600nm','650nm','700nm','750nm','800nm','850nm')
title('(a)N=4');
ylabel(' Transmission');
xlabel('¦Ë[nm]');
axis([1250 1750,-0.2 1.2])
subplot(1,2,2)
plot(g,t_1,'k',g,t_2,'r',g,t_3,'g',g,t_4,'b',g,t_5,'y',g,t_6,'m')
legend('600nm','650nm','700nm','750nm','800nm','850nm')
title('(b)N=5');
ylabel(' Transmission');
xlabel('¦Ë[nm]');
axis([1250 1750,-0.2 1.2])


