clear all
clc
cd /Users/Daniel/Downloads/ACMsimulation/NaCl/naclp6
tt=1; %time step 1fs;
tf=1000; %output frequency;
filename = 'Rg_nvt8.data';%change data file name;

delimiterIn = ' ';
A = importdata(filename);
B = A.data;

xp=B(1,2); %Number of molecules lines;
l=size(B,1); %length of the matrix;

t=(l/(xp+1))*tt*tf; % time of the simulaiton in fs;
s=t/tf; % 1000 polymer chains or chunks;

        
for o=1:s
    i=(xp+1)*(o-1);
    for j=1:xp
        k=i+j+1;
        C(j,o)=B(k,2);
    end
end

l=size(C,2);
for i=1:l
    T(1,i)=(i*tt*tf)/1000000; %time in ns;
end

figure;
plot(T(1,:),C(93,:),'r');
title('PAH200 Rg vs Time');
xlabel('time (ns)');
ylabel('PAH200 Rg (Å)');

result(:,1)=transpose(T(1,1:l));
result(:,2)=transpose(C(93,1:l));
