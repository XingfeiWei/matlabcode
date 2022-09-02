clear all;
clc
cd ~/Desktop/
filename = 'rdf_MC_d10_RDF.txt';%change data file name;
delimiterIn = ' ';
A = importdata(filename);
B=A.data;
ld=B(1,2);
nd=ld*2+1;
n=size(B,1)/nd;
for i=1:n
    for j = 1:ld
        k=(i-1)*nd+j*2+1;
        C(j,i)=B(k,1);
    end
end

for i=1:1
    for j = 1:ld
        k=(i-1)*nd+j*2;
        R(j,i)=B(k,2);
    end
end
%%
clear RDF
RDF(:,1) = C(:,1);
for i = 1 :ld
RDF(i,2) = mean(C(i,n-80:n-60));
end
for i = 1 :ld
RDF(i,3) = mean(C(i,n-60:n-40));
end
for i = 1 :ld
RDF(i,4) = mean(C(i,n-40:n-20));
end
for i = 1 :ld
RDF(i,5) = mean(C(i,n-20:n));
end

figure; plot(R,RDF);
