clear all
clc
cd Rg1
tt=0.0005; %time step 1fs;
tf=10000; %output frequency;
rgfile1     =        'Rg_A500Rg1.txt';%change data file name;
logfile1=importdata('log.A500Rg1_1',' ',1370);
delimiterIn = ' ';
A = importdata(rgfile1);
B = A.data;
l=size(B,1); %length of the matrix;
for i=1:l/2
    k=i*2;
    C(i,1)=B(k,2)*8.518; %convert to A
    t(i,1)=i*tt*tf*3.03/1000; %convert to ns
end

figure;
subplot(2,1,1);
plot(t,C);
title('Rg vs Time');
xlabel('time (ns)');
ylabel('Rg (Å)');
hold on

ac=logfile1.data;
tlog=ac(:,1)*tt*3.03/1000; %convert to ns
PE=ac(:,3);

subplot(2,1,2);
plot(tlog,PE,'r');
title('PotEnergy vs Time');
xlabel('Time (ns)');
ylabel('PetEnergy (kcal/mol)');
%%
result(1,1)=mean(C(l/2-5000:l/2,1))

