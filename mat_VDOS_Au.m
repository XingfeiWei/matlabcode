%% VDOS.m

clear all;

filename = 'av1.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data1 finished!');

dt=D(2,1)-D(1,1); %step size
timestep=0.25*dt; % in fs
tmax=0.5;% time length in ns

nmax=tmax/timestep*1e6+1;% maximum number of data
l=size(D,1);

fileID = fopen('timea.txt','w');
fprintf(fileID,'%6s \n','#time in ps');
for n=1:nmax
t(n,1)=timestep*(n-1)/1000;% in ps
fprintf(fileID,'%6.8f \n',t(n,1));
end
fclose(fileID);
fprintf('output time finished!');

vx=D(:,2);
vy=D(:,3);
vz=D(:,4);

parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv1(n,:)=[tem,temx,temy,temz];
end
%%
fileID1 = fopen('vv1a.txt','w');
fprintf(fileID1,'%s \n','#VACF of atom1 vv xx yy zz');
for i=1:nmax
fprintf(fileID1,'%e ',vv1(i,1));
fprintf(fileID1,'%e ',vv1(i,2));
fprintf(fileID1,'%e ',vv1(i,3));
fprintf(fileID1,'%e \n',vv1(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom1 finished!\n');

%% data2
filename = 'av2.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data2 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv2(n,:)=[tem,temx,temy,temz];
end
fileID1 = fopen('vv2a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom2');
for i=1:nmax
fprintf(fileID1,'%e ',vv2(i,1));
fprintf(fileID1,'%e ',vv2(i,2));
fprintf(fileID1,'%e ',vv2(i,3));
fprintf(fileID1,'%e \n',vv2(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom2 finished!');

%data3
filename = 'av3.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data3 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv3(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv3a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom3');
for i=1:nmax
fprintf(fileID1,'%e ',vv3(i,1));
fprintf(fileID1,'%e ',vv3(i,2));
fprintf(fileID1,'%e ',vv3(i,3));
fprintf(fileID1,'%e \n',vv3(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom3 finished!');

%data4
filename = 'av4.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data4 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv4(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv4a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom4');
for i=1:nmax
fprintf(fileID1,'%e ',vv4(i,1));
fprintf(fileID1,'%e ',vv4(i,2));
fprintf(fileID1,'%e ',vv4(i,3));
fprintf(fileID1,'%e \n',vv4(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom4 finished!');

%data5
filename = 'av5.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data5 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv5(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv5a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom5');
for i=1:nmax
fprintf(fileID1,'%e ',vv5(i,1));
fprintf(fileID1,'%e ',vv5(i,2));
fprintf(fileID1,'%e ',vv5(i,3));
fprintf(fileID1,'%e \n',vv5(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom5 finished!');

%data6
filename = 'av6.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data6 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv6(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv6a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom6');
for i=1:nmax
fprintf(fileID1,'%e ',vv6(i,1));
fprintf(fileID1,'%e ',vv6(i,2));
fprintf(fileID1,'%e ',vv6(i,3));
fprintf(fileID1,'%e \n',vv6(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom6 finished!');

%data7
filename = 'av7.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data7 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv7(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv7a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom7');
for i=1:nmax
fprintf(fileID1,'%e ',vv7(i,1));
fprintf(fileID1,'%e ',vv7(i,2));
fprintf(fileID1,'%e ',vv7(i,3));
fprintf(fileID1,'%e \n',vv7(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom7 finished!');

%data8
filename = 'av8.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data8 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv8(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv8a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom8');
for i=1:nmax
fprintf(fileID1,'%e ',vv8(i,1));
fprintf(fileID1,'%e ',vv8(i,2));
fprintf(fileID1,'%e ',vv8(i,3));
fprintf(fileID1,'%e \n',vv8(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom8 finished!');

%data9
filename = 'av9.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data9 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv9(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv9a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom9');
for i=1:nmax
fprintf(fileID1,'%e ',vv9(i,1));
fprintf(fileID1,'%e ',vv9(i,2));
fprintf(fileID1,'%e ',vv9(i,3));
fprintf(fileID1,'%e \n',vv9(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom9 finished!');

%data10
filename = 'av10.txt';%change data file name;
delimiterIn = ' ';
headerlinesIn = 2;
C = importdata(filename,delimiterIn,headerlinesIn);
D = C.data;
fprintf('load data10 finished!');
vx=D(:,2);
vy=D(:,3);
vz=D(:,4);
parfor (n=1:nmax, 12)
tem=0;temx=0;temy=0;temz=0;
for i=1:l-n+1
	tem = vx(i,1)*vx(i+n-1,1)+vy(i,1)*vy(i+n-1,1)+vz(i,1)*vz(i+n-1,1)+tem;
    temx = vx(i,1)*vx(i+n-1,1)+temx;
    temy = vy(i,1)*vy(i+n-1,1)+temy;
    temz = vz(i,1)*vz(i+n-1,1)+temz;
end
tem=tem/(l-n+1);temx=temx/(l-n+1);temy=temy/(l-n+1);temz=temz/(l-n+1);
vv10(n,:)=[tem,temx,temy,temz];
end

fileID1 = fopen('vv10a.txt','w');
fprintf(fileID1,'%6s \n','#VACF of atom10');
for i=1:nmax
fprintf(fileID1,'%e ',vv10(i,1));
fprintf(fileID1,'%e ',vv10(i,2));
fprintf(fileID1,'%e ',vv10(i,3));
fprintf(fileID1,'%e \n',vv10(i,4));
end
fclose(fileID1);
fprintf('output VACF of atom10 finished!');

%calcualte the mean VACF
fileID1 = fopen('vva.txt','w');
fprintf(fileID1,'%6s \n','#VACF of mean');
for i=1:nmax
vv(i,1)=(vv1(i,1)+vv2(i,1)+vv3(i,1)+vv4(i,1)+vv5(i,1)+vv6(i,1)+vv7(i,1)+vv8(i,1)+vv9(i,1)+vv10(i,1))/10;
vv(i,2)=(vv1(i,2)+vv2(i,2)+vv3(i,2)+vv4(i,2)+vv5(i,2)+vv6(i,2)+vv7(i,2)+vv8(i,2)+vv9(i,2)+vv10(i,2))/10;
vv(i,3)=(vv1(i,3)+vv2(i,3)+vv3(i,3)+vv4(i,3)+vv5(i,3)+vv6(i,3)+vv7(i,3)+vv8(i,3)+vv9(i,3)+vv10(i,3))/10;
vv(i,4)=(vv1(i,4)+vv2(i,4)+vv3(i,4)+vv4(i,4)+vv5(i,4)+vv6(i,4)+vv7(i,4)+vv8(i,4)+vv9(i,4)+vv10(i,4))/10;
fprintf(fileID1,'%e ',vv(i,1));
fprintf(fileID1,'%e ',vv(i,2));
fprintf(fileID1,'%e ',vv(i,3));
fprintf(fileID1,'%e \n',vv(i,4));
end
fclose(fileID1);
fprintf('output VACF of all atoms average finished!');

%Fourier transform
Fs=1/(timestep*1e-15);%frequency per sample fs
L=length(vv);

nf=2^nextpow2(L);
Fvv=fft(vv(:,1),nf,1);
Fvvx=fft(vv(:,2),nf,1);
Fvvy=fft(vv(:,3),nf,1);
Fvvz=fft(vv(:,4),nf,1);

for i=1:nf/2+1 
	Fre(i,1)=Fs*(i-1)/nf/1e12; %to THz
    VDOS(i,1)=abs(Fvv(i,1));
    VDOS(i,2)=abs(Fvvx(i,1));
    VDOS(i,3)=abs(Fvvy(i,1));
    VDOS(i,4)=abs(Fvvz(i,1));
end 

fileID = fopen('VDOSa.txt','w');
	fprintf(fileID,'Frequency VDOS \n');
for i=1:nf/2+1 
	fprintf(fileID,'%e ',Fre(i,1));	
    fprintf(fileID,'%e ',VDOS(i,1));
	fprintf(fileID,'%e ',VDOS(i,2));
	fprintf(fileID,'%e ',VDOS(i,3));
    fprintf(fileID,'%e \n',VDOS(i,4));
end
fclose(fileID);
fprintf('output VDOS finished!');

h=figure;
plot(Fre,VDOS);title('VDOS of Au');xlabel('Frequency (THz)');ylabel('VDOS');
ylim([0,5e-1]);xlim([0,250]);
savefig(h,'VDOS_Au.fig');

Fvv1=fft(vv1(:,1),nf,1);
Fvv2=fft(vv2(:,1),nf,1);

for i=1:nf/2+1 
    VDOS1(i,1)=abs(Fvv1(i,1));
    VDOS2(i,1)=abs(Fvv2(i,1));
end 
h1=figure;
plot(Fre,VDOS1,Fre,VDOS2);title('VDOS of Au1 Au2');xlabel('Frequency (THz)');ylabel('VDOS');
ylim([0,5e-1]);xlim([0,250]);
savefig(h1,'VDOS_Au1_Au2.fig');


