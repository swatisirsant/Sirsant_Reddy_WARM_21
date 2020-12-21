function lcc = TL_lcc_gen_exp(Dia,T,length)
ch = [0 25.4 50.8 76.2 101.6 152.4 203.2 254 304.8 355.6 406.4 457.2 508 558.8 609.6];    
cost=  [0 2.1	5.2	8.9	12.01	17.35	24.93	34.12	52.51	63.22	93.8	134.5	175.8	306.12	556.44]; % pipe+installation unit cost
No=[0 1.3 1.05 0.81 0.58 0.41 0.25 0.15 0.1 0.08 0.06 0.05 0.04 0.03 0.02]; %initial break rate
br_cost=[0 56.68	61.84	64.41	83.74	135.26	199.68	270.5	360.7	463.77	579.71	734.3	966.18	1127.21	1288.24]; % break cost per break
rep_cost=[0 6.7 8.12 9.8 12.02 17.76 25.67 34.18 43.98 55.1 67.4 81.9	103.78 119.88 136]; % replacement cost per m
A=0.08; %break growth rate coefficient
y=50; %remaining planning horizon in years
netcst=0;
cst =zeros(1,16);
br_int=zeros(1,16);
b_cst=zeros(1,16);
r_cst=zeros(1,16);
for j=1:16
    for k=1:15
       if round(Dia(j)*1000)==round(ch(k)*1000);
            cst(j)=length(j)*cost(k);
            br_int(j)=No(k);
            b_cst(j)=br_cost(k);
            r_cst(j)=rep_cost(k);
       end
    end
end
for i=1:16
    netcst=netcst+cst(i)*(1.04)^(T(i)-10)*1.08^(y-(T(i)-10)); %cost of pipe +installation
end


CR=0.0048*Dia+6.6805; %cost ratio
BRth=zeros(1,16);
for i=1:16
    BRth(i)=log(1.04)/log(1+CR(i)/length(i));
end
for i=1:16
    if Dia(i)==0
        s(i)=0;
        k(i)=0;
    else
        s(i)=fix(log(BRth(i)/br_int(i))/A); % s is service life of pipe i and A is break growth rate coefficient
        k(i)=fix((y-(T(i)-10))/s(i)); % y is planning horizon in years, k is the number of break repairs needed
    end
end
BR=zeros(1,16);
RC=zeros(1,16);
AR=zeros(1,16);
for i=1:16
    if k(i)>0
        for m=1:k(i)
            for r=1:s(i)
                t=((m-1)*s(i))+r;
                BR(i)=BR(i)+(br_int(i)*exp(A*t)*b_cst(i)/1.04^(t+T(i)-10));
            end
            RC(i)=RC(i)+(r_cst(i)*length(i)/1.04^(m*s(i)+T(i)-10));
        end
    end
    if s(i)>0
        for r=(k(i)*s(i)+1):(y-(T(i)-10))
            AR(i)=AR(i)+(br_int(i)*exp(A*(r-k(i)*s(i)))*b_cst(i)/1.04^(r+T(i)-10));
        end
    end
end
BR_t=0;
RC_t=0;
AR_t=0;
for i=1:16
    BR_t=BR(i)*(1.08)^(y-(T(i)-10));
    RC_t=RC(i)*(1.08)^(y-(T(i)-10));
    AR_t=AR(i)*(1.08)^(y-(T(i)-10));
end
lcc=netcst+BR_t+RC_t+AR_t;