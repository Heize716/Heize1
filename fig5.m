%QÖµ
%×¢ÒâÐè·Ö°å¿é¼ÆËã
%Note that calculations need to be divided into sections
%% (BA)^{4}(BCB) ^{1}(AB)^{4}ÐÍ
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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

%ÖÐÐÄ²¨³¤
data=[g;t1;t2;t3;t4;t5;t6]'
[M1,I1]=max(t1)     
t1(I1)              %ÖÐÐÄ²¨³¤·åÖµ
P1=data(I1,1);
W1=P1               %ÖÐÐÄ²¨³¤

[M2,I2]=max(t2)
t2(I2)
P2=data(I2,1);
W2=P2

[M3,I3]=max(t3)
t3(I3)
P3=data(I3,1);
W3=P3

[M4,I4]=max(t4)
t4(I4)
P4=data(I4,1);
W4=P4

[M5,I5]=max(t5)
t5(I5)
P5=data(I5,1);
W5=P5

[M6,I6]=max(t6)
t6(I6)
P6=data(I6,1);
W6=P6


%°ë·å¿í¶È
[A1, A_1]=min(abs(t1(1:I1)-0.5))
[A2, A_2]=min(abs(t1(A_1+1:end)-0.5))
x1=t1(A_1)              %°ë·åÖµ(×ó)
y1=t1(A_1+A_2)          %°ë·åÖµ(ÓÒ)
X1=data(A_1,1);         %°ë·åÖµ×ó±ß²¨³¤
Y1=data(A_1+A_2,1);     %°ë·åÖµ×ó±ß²¨³¤

[B1, B_1]=min(abs(t2(1:I2)-0.5))
[B2, B_2]=min(abs(t2(B_1+1:end)-0.5))
x2=t2(B_1)
y2=t2(B_1+B_2)
X2=data(B_1,1);
Y2=data(B_1+B_2,1);

[C1, C_1]=min(abs(t3(1:I3)-0.5))
[C2, C_2]=min(abs(t3(C_1+1:end)-0.5))
x3=t3(C_1)
y3=t3(C_1+C_2)
X3=data(C_1,1);
Y3=data(C_1+C_2,1);

[D1, D_1]=min(abs(t4(1:I4)-0.5))
[D2, D_2]=min(abs(t4(D_1+1:end)-0.5))
x4=t4(D_1)
y4=t4(D_1+D_2)
X4=data(D_1,1);
Y4=data(D_1+D_2,1);

[E1, E_1]=min(abs(t5(1:I5)-0.5))
[E2, E_2]=min(abs(t5(E_1+1:end)-0.5))
x5=t5(E_1)
y5=t5(E_1+E_2)
X5=data(E_1,1);
Y5=data(E_1+E_2,1);

[F1, F_1]=min(abs(t6(1:I6)-0.5))
[F2, F_2]=min(abs(t6(F_1+1:end)-0.5))
x6=t6(F_1)
y6=t6(F_1+F_2)
X6=data(F_1,1);
Y6=data(F_1+F_2,1);


%QÖµ
Q1=W1/(Y1-X1);
Q2=W2/(Y2-X2);
Q3=W3/(Y3-X3);
Q4=W4/(Y4-X4);
Q5=W5/(Y5-X5);
Q6=W6/(Y6-X6);

%% (BA)^{5}(BCB) ^{1}(AB)^{5}ÐÍ
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

t1=[];
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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
g=linspace(1250,1750,10000);%nm
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

%ÖÐÐÄ²¨³¤
data=[g;t1;t2;t3;t4;t5;t6]'
[M1,I1]=max(t1)     
t1(I1)              %ÖÐÐÄ²¨³¤·åÖµ
P1=data(I1,1);
W1=P1               %ÖÐÐÄ²¨³¤

[M2,I2]=max(t2)
t2(I2)
P2=data(I2,1);
W2=P2

[M3,I3]=max(t3)
t3(I3)
P3=data(I3,1);
W3=P3

[M4,I4]=max(t4)
t4(I4)
P4=data(I4,1);
W4=P4

[M5,I5]=max(t5)
t5(I5)
P5=data(I5,1);
W5=P5

[M6,I6]=max(t6)
t6(I6)
P6=data(I6,1);
W6=P6


%°ë·å¿í¶È
[A1, A_1]=min(abs(t1(1:I1)-0.5))
[A2, A_2]=min(abs(t1(A_1+1:end)-0.5))
x1=t1(A_1)              %°ë·åÖµ(×ó)
y1=t1(A_1+A_2)          %°ë·åÖµ(ÓÒ)
X1=data(A_1,1);         %°ë·åÖµ×ó±ß²¨³¤
Y1=data(A_1+A_2,1);     %°ë·åÖµ×ó±ß²¨³¤

[B1, B_1]=min(abs(t2(1:I2)-0.5))
[B2, B_2]=min(abs(t2(B_1+1:end)-0.5))
x2=t2(B_1)
y2=t2(B_1+B_2)
X2=data(B_1,1);
Y2=data(B_1+B_2,1);

[C1, C_1]=min(abs(t3(1:I3)-0.5))
[C2, C_2]=min(abs(t3(C_1+1:end)-0.5))
x3=t3(C_1)
y3=t3(C_1+C_2)
X3=data(C_1,1);
Y3=data(C_1+C_2,1);

[D1, D_1]=min(abs(t4(1:I4)-0.5))
[D2, D_2]=min(abs(t4(D_1+1:end)-0.5))
x4=t4(D_1)
y4=t4(D_1+D_2)
X4=data(D_1,1);
Y4=data(D_1+D_2,1);

[E1, E_1]=min(abs(t5(1:I5)-0.5))
[E2, E_2]=min(abs(t5(E_1+1:end)-0.5))
x5=t5(E_1)
y5=t5(E_1+E_2)
X5=data(E_1,1);
Y5=data(E_1+E_2,1);

[F1, F_1]=min(abs(t6(1:I6)-0.5))
[F2, F_2]=min(abs(t6(F_1+1:end)-0.5))
x6=t6(F_1)
y6=t6(F_1+F_2)
X6=data(F_1,1);
Y6=data(F_1+F_2,1);


%QÖµ
QQ1=W1/(Y1-X1);
QQ2=W2/(Y2-X2);
QQ3=W3/(Y3-X3);
QQ4=W4/(Y4-X4);
QQ5=W5/(Y5-X5);
QQ6=W6/(Y6-X6);


%% QÖµ
X=[600 650 700 750 800 850 ]
Y=[Q1 Q2 Q3 Q4 Q5 Q6]
Z=[QQ1 QQ2 QQ3 QQ4 QQ5 QQ6]
figure
plot(X,Y,'-o',X,Z,'-*')
title('N-Q');
ylabel('Q');
xlabel('¦Ë[nm]');
legend('N=4','N=5')