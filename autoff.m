
clear all;
fname='FF';
delimiterIn = ' ';
headerlinesIn = 2;
a=importdata(fname,delimiterIn,headerlinesIn);
b=a.data;
l=size(b,1);
for i=1:l
    P(i,1)=b(i,2); % potential in kcal/mol
    F(i,1)=b(i,5); % Force in z direction
end

nmax=(l-1)/2+1;
%nmax=1000;
timestep=1;

dt=1;
delt=dt*timestep/1e15;%in s
kB=1.38064852e-23; % in J/K
T=300;% K
unit=4186.6/6.022140857e23/1e-10;%kcal/mol/A to N
R=5.427979; %Radius in A
L=27.016529; %Length in A
S=3.1415926*2*R*L*1e-20; %m^2

fileID = fopen('ff0.txt','w');
fprintf(fileID,'%6s %12s %12s\n','#time','<F(0).F(t)>','Friction coeff in N.s/m^3, normal+Trap+Simp');

FF=zeros(nmax,1);
CF0=zeros(nmax,1);
CF1=zeros(nmax,1);
CF2=zeros(nmax,1);
temG=zeros(nmax,1);

parfor (n=1:nmax, 12); %% only this for loop is changed to a parallel for loop
tem=0;
for i=1:l-n+1
tem = F(i,1)*F(i+n-1,1)+tem;
end
tem = tem/(l-n+1);
FF(n,1)=tem;
end

for i=1:nmax
    t(i,1)=(i-1)/1000; % time in ps
if i==1
CF0(i,1)=FF(i,1)*unit*unit*delt;
CF1=0;
CF2=0;
else
if i==2
    CF0(i,1)=CF0(i-1,1)+FF(i,1)*unit*unit*delt; % One point integration
	CF1(i,1)=CF1(i-1,1)+0.5*(FF(i-1,1)+FF(i,1))*unit*unit*delt; %%Trapzoidal rule
	CF2(i,1)=CF2(i-1,1)+0.5*(FF(i-1,1)+FF(i,1))*unit*unit*delt; %%Simpson rule
else
CF0(i,1)=CF0(i-1,1)+FF(i,1)*unit*unit*delt; % One point integration
CF1(i,1)=CF1(i-1,1)+0.5*(FF(i-1,1)+FF(i,1))*unit*unit*delt; %%Trapzoidal rule
CF2(i,1)=CF2(i-1,1)+(FF(i-2,1)+4*FF(i-1,1)+FF(i,1))*unit*unit*delt/6; %%Simpson rule
end
end
coef0(i,1)=CF0(i,1)/kB/T/S;
coef1(i,1)=CF1(i,1)/kB/T/S;
coef2(i,1)=CF2(i,1)/kB/T/S;
fprintf(fileID,'%.8f %12.8f %12.8f %12.8f %12.8f\n',t(i,1),FF(i,1),coef0(i,1),coef1(i,1),coef2(i,1));
end
fclose(fileID);
