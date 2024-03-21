% (BA)^{4}(BCB) ^{1}(AB)^{4}ÐÍ
clear all;k=4;n=1;
n0=1;na=3.230;nb=1.35;nc=1.0;
c1=0/6; %½Ç¶È¿Éµ÷
ca=asin(n0*sin(c1)/na);cb=asin(na*sin(ca)/nb);cc=asin(nb*sin(cb)/nc);
c6=c1;

d0=87.5;da=1480/(4*na);db=1480/(4*nb);dc1=650;dc2=675;dc3=700;dc4=725;dc5=750;dc6=775;d4=87.5; %µ¥Î»nm
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
g=linspace(1250,1700,10000);%nm
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
g=linspace(1250,1700,10000);%nm
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
g=linspace(1250,1700,10000);%nm
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
g=linspace(1250,1700,10000);%nm
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
g=linspace(1250,1700,10000);%nm
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
g=linspace(1250,1700,10000);%nm
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


%R-ÆµÂÊ
cc=3e17
f=cc./g*10^(-12)
figure
subplot(1,2,1)
plot(f,1-t1,'k',f,1-t2,'r',f,1-t3,'g',f,1-t4,'b',f,1-t5,'y',f,1-t6,'m')
legend('600nm','650nm','700nm','750nm','800nm','850nm')
title('(a)Resonance peak reflection frequency shift diagram');
ylabel(' Reflectivity');
xlabel('Frequency[THz]');
axis([1.9e2 2.2e2,-0.2 1.2])

%ÖÐÐÄ²¨³¤
data=[g;t1;t2;t3;t4;t5;t6]'
[M1,I1]=max(t1)
P1=data(I1,1);

[M2,I2]=max(t2)
P2=data(I2,1);

[M3,I3]=max(t3)
P3=data(I3,1);

[M4,I4]=max(t4)
P4=data(I4,1);

[M5,I5]=max(t5)
P5=data(I5,1);

[M6,I6]=max(t6)
P6=data(I6,1);

%²¨³¤-Ç»³¤
cc=3e17
X=[600 650 700 750 800 850]
Y=[cc./P1*10^(-12) cc./P2*10^(-12) cc./P3*10^(-12) cc./P4*10^(-12) cc./P5*10^(-12) cc./P6*10^(-12)];
t=polyfit(Y,X,1);
subplot(1,2,2)
plot(Y,X,'o',Y,polyval(t,Y))
%legend('Simulation','Linear fit','600nm','650nm''700nm','750nm','800nm','850nm','850nm')
legend('Simulation','Linear fit')
title('(b)Fitting curve');
ylabel(' Cavity Length [nm]');
xlabel('Frequency[THz]');
axis([1.98e2 2.12e2 ,580 880])
