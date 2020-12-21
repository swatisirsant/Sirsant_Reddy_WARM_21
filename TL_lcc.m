function lc_cost = TL_lcc(dia1,length)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% wdsfile='TL.inp';
% addpath('G:\matlab_codes\Epanet files');
% epanetloadfile(wdsfile);
ch = [25.4 50.8 76.2 101.6 152.4 203.2 254 304.8 355.6 406.4 457.2 508 558.8 609.6];    
cost=  [2.1	5.2	8.9	10.84	16.15	23.33	31.07	39.99	50.09	61.27	74.45	94.35	108.99	123.63]; % pipe+installation unit cost
No=[1.3 1.05 0.81 0.58 0.41 0.25 0.15 0.1 0.08 0.06 0.05 0.04 0.03 0.02]; %initial break rate
br_cost=[56.68	61.84	64.41	83.74	135.26	199.68	270.5	360.7	463.77	579.71	734.3	966.18	1127.21	1288.24]; % break cost per break
rep_cost=[6.7 8.12 9.8 12.02 17.76 25.67 34.18 43.98 55.1 67.4 81.9	103.78 119.88 136]; % replacement cost per m
A=0.08; %break growtn rate coefficient
y=50; %planning horizon in years
cst =zeros(1,8);
br_int=zeros(1,8);
b_cst=zeros(1,8);
r_cst=zeros(1,8);
for j=1:8
    for k=1:14
       if round(dia1(j)*1000)==round(ch(k)*1000);
            cst(j)=length(j)*cost(k);
            br_int(j)=No(k);
            b_cst(j)=br_cost(k);
            r_cst(j)=rep_cost(k);
       end
    end
end
netcst=sum(cst); %cost of pipe +installation


CR=0.0048*dia1+6.6805; %cost ratio
BRth=zeros(1,8);
for i=1:8
    BRth(i)=log(1.04)/log(1+CR(i)/length(i));
end
for i=1:8
    s(i)=fix(log(BRth(i)/br_int(i))/A); % s is service life of pipe i and A is break growth rate coefficient
    k(i)=fix(y/s(i)); % y is planning horizon in years, k is the number of break repairs needed
end
BR=zeros(1,8);
RC=zeros(1,8);
AR=zeros(1,8);
for i=1:8
    if k(i)>0
        for m=1:k(i)
            for r=1:s(i)
                t=((m-1)*s(i))+r;
                BR(i)=BR(i)+(br_int(i)*exp(A*t)*b_cst(i)/1.04^t);
            end
            RC(i)=RC(i)+(r_cst(i)*length(i)/1.04^(m*s(i)));
        end
    end
    for r=(k(i)*s(i)+1):y
        AR(i)=AR(i)+(br_int(i)*exp(A*(r-k(i)*s(i)))*b_cst(i)/1.04^r);
    end
end

BR_t=sum(BR);
RC_t=sum(RC);
AR_t=sum(AR);
lc_cost=(netcst+BR_t+RC_t+AR_t)*1.08^y;