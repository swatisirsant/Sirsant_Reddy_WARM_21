clear all;
clc;
wdsfile='TL_exp1.inp';
addpath('D:\Epanet_codes');
epanetloadfile(wdsfile);
Npipes = 16;
Nnodes = getdata('EN_NODECOUNT')-1;
Length = getdata('EN_LENGTH');
Demand = getdata('EN_BASEDEMAND');
ele = getdata('EN_ELEVATION');
Hmin=ones(1,Nnodes)*30;
Dia=[609.6	355.6	457.2	558.8	152.4	457.2	254	609.6];
% Dia1=[609.6	0	609.6	406.4	203.2	0	0	0	609.6	0	609.6	609.6	152.4	508	0	609.6];
dia=[25.4 50.8 76.2 101.6 152.4 203.2 254 304.8 355.6 406.4 457.2 508 558.8 609.6];
I=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
Demand1=Demand*1.03^20;
setdata('EN_BASEDEMAND',Demand1);
N_poss=16; %max possible no of pipes that can be added
for N=1:N_poss
    combos = randsample(dia,N*N*10,true);
    combos=reshape(combos,N*10,N);
    [s1 s2]=size(combos);
    choices=randsample(I,N*N*10,true);
    choices=reshape(choices,N*10,N); %  pipes selected for parallelising, N at a time
    [p1,p2]=size(choices);
    epanetclose();
    epanetloadfile(wdsfile);
    for i=1:p1
        cost1=ones(1,s1)*10^40;
        res1=zeros(1,s1);
        d=zeros(s1,16);
        k=1;
        T(i,:)=zeros(1,16);
        for j=1:16
            if k>N
                break;
            end
            if j==choices(i,k)
                for l=1:s1
                    d(l,j)=combos(l,k);
                end
                k=k+1;
                T(i,j)=20;
            end
        end
        for j=1:s1
        D=d(j,:);
        DD=[Dia D];
        
        setdata('EN_DIAMETER',DD);
        setdata('EN_BASEDEMAND',Demand1);
        ENsolveH();
        h=getdata('EN_PRESSURE');
%         pause
        flow=getdata('EN_FLOW');
        for l=1:Nnodes
            if h(l)<Hmin(l)
                DH(l)=Hmin(l)-h(l);
            else
                DH(l)=0;
            end
        end
        DHmax=max(DH);
        Penalty=DHmax*10^20;
        tot_h=h+ele;
        cost1(j)= TL_lcc_gen_exp(D,T(i,:),Length)+Penalty;
        res1(j)=Res_TL_exp2(tot_h,Demand1,flow);
        if res1(j)<0.8
            Penaltyr=(0.8-res1(j))*10^20;
        else
            Penaltyr=0;
        end
        cost1(j)=cost1(j)+Penaltyr;
        end
        [cmin1(i),A1]=min(cost1(1,:));
        res_m1(i)=res1(A1);
        D_best1(i,:)=d(A1,:);
        T_best1(i,:)=T(i,:);
    end
    [cmin1_1(N),A1]=min(cmin1);
    res_m1_1(N)=res_m1(A1);
    D_best1_1(N,:)=D_best1(A1,:);
    T_best1_1(N,:)=T_best1(A1,:);
    clear cost1;
    clear res1;
    clear D
    clear T;
end

%%stage2
Demand2=Demand*1.03^30;
epanetclose();
epanetloadfile(wdsfile);
% setdata('EN_BASEDEMAND',Demand2);
for N=1:N_poss 
    for n=1:N-1 %start at n=1, since at stage 1 minimum 1 pipe need to be added since without adding any pipe at stage 1, min head requirement is not satisfied
        dd1=D_best1_1(n,:);
        m=(dd1~=0);
        I1=[];
        for q=1:16
            for r=1:16
                if m(q)==0 && q==I(r)
                    I1=[I1 I(r)];
                end
            end
        end
%         [s,s1]=size(I1);
%         p1=factorial(s1)/(factorial(N-n)*factorial(s1-N+n));
        choices=randsample(I1,(N-n)*(N-n)*10,true);
        choices=reshape(choices,10*(N-n),(N-n));
        [p1,p2]=size(choices);
%         choices=combntns(I1,N-n);
%         combos = dec2base((0:14^(N-n)-1)', 14) - '0';
%         combos = sort(combos, 2);
%         combos = unique(combos, 'rows');
        combos=randsample(dia,(N-n)*(N-n)*10,true);
        combos=reshape(combos,10*(N-n),(N-n));
        [s1,s2]=size(combos);
%         for i=1:s1
%             for j=1:s2
%                 if combos(i,j)>9
%                     combos(i,j)=combos(i,j)-7;
%                 end
%             end
%         end
%         combos = dia(1 + combos);
%         if (N-n)==1
%             combos=combos';
%         end
        d=zeros(s1,16);
        for i=1:s1
            d(i,:)=dd1;
        end
        for i=1:p1
%                 d=dd1;
                cost2=ones(1,s1)*10^40;
                res2_1=zeros(1,s1);
                res2_2=zeros(1,s1);
                k=1;
                T(i,:)=m*20; % when m=1, T is 20 yrs
                T1(i,:)=T(i,:);
                for j=1:16
                    if k>(N-n)
                        break;
                    end
                    if j==choices(i,k)
                        for l=1:s1
                            d(l,j)=combos(l,k);
                        end
                        k=k+1;
                        T(i,j)=30;
                    end
                end
                for j=1:s1
                D=d(j,:);
                DD=[Dia dd1];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand1);
                ENsolveH();
                h1=getdata('EN_PRESSURE');
        %         pause
                flow1=getdata('EN_FLOW');
                DD=[Dia D];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand2);
                ENsolveH();
                h2=getdata('EN_PRESSURE');
                flow2=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h1(l)<Hmin(l)
                        DH1(l)=Hmin(l)-h1(l);
                    else
                        DH1(l)=0;
                    end
                end
                DHmax1=max(DH1);
                Penalty1=DHmax1*10^20;
                tot_h1=h1+ele;
                for l=1:Nnodes
                    if h2(l)<Hmin(l)
                        DH2(l)=Hmin(l)-h2(l);
                    else
                        DH2(l)=0;
                    end
                end
                DHmax2=max(DH2);
                Penalty2=DHmax2*10^20;
                tot_h2=h2+ele;
                cost2(j)= TL_lcc_gen_exp(D,T(i,:),Length)+Penalty1+Penalty2;
                res2_1(j)=Res_TL_exp2(tot_h1,Demand1,flow1);
                res2_2(j)=Res_TL_exp2(tot_h2,Demand2,flow2);
                if res2_1(j)<0.8
                    Penaltyr1=(0.8-res2_1(j))*10^20;
                else
                    Penaltyr1=0;
                end
                if res2_2(j)<0.8
                    Penaltyr2=(0.8-res2_2(j))*10^20;
                else
                    Penaltyr2=0;
                end
                cost2(j)=cost2(j)+Penaltyr1+Penaltyr2;
                end
                [cmin2(i),A1]=min(cost2(1,:));
                res_m2(i)=res2_2(A1);
                D_best2(i,:)=d(A1,:);
                T_best2(i,:)=T(i,:);
        end
            [cmin2_1(N),A1]=min(cmin2);
            res_m2_1(N)=res_m2(A1);
            D_best2_1(N,:)=D_best2(A1,:);
            T_best2_1(N,:)=T_best2(A1,:);
            clear cost2;
            clear res2;
            clear D;
            clear T;
    end
            D_best2_1(N,:)=D_best1_1(N,:);
            T_best2_1(N,:)=T_best1_1(N,:);
            DD=[Dia D_best2_1(N,:)];
%             epanetclose();
%             epanetloadfile(wdsfile);
            setdata('EN_BASEDEMAND',Demand1);
            setdata('EN_DIAMETER',DD);
            ENsolveH();
            h1=getdata('EN_PRESSURE');
            flow1=getdata('EN_FLOW');
%             epanetclose();
%             epanetloadfile(wdsfile);
            setdata('EN_DIAMETER',DD);
            setdata('EN_BASEDEMAND',Demand2);
            ENsolveH();
            h2=getdata('EN_PRESSURE');
            flow2=getdata('EN_FLOW');
            for l=1:Nnodes
                if h1(l)<Hmin(l)
                    DH1(l)=Hmin(l)-h1(l);
                else
                    DH1(l)=0;
                end
            end
            DHmax1=max(DH1);
            Penalty1=DHmax1*10^20;
            tot_h1=h1+ele;
            for l=1:Nnodes
                if h2(l)<Hmin(l)
                    DH2(l)=Hmin(l)-h2(l);
                else
                    DH2(l)=0;
                end
            end
            DHmax2=max(DH2);
            Penalty2=DHmax2*10^20;
            tot_h2=h2+ele;
            cmin2_1(N)=TL_lcc_gen_exp(D_best2_1(N,:),T_best2(N,:),Length)+Penalty1+Penalty2;
            res2=Res_TL_exp2(tot_h1,Demand1,flow1);
            if res2<0.8
                Penaltyr1=(0.8-res2)*10^20;
            else
                Penaltyr1=0;
            end
            res_m2(N)=Res_TL_exp2(tot_h2,Demand2,flow2);
            if res_m2(N)<0.8
                Penaltyr2=(0.8-res_m2(N))*10^20;
            else
                Penaltyr2=0;
            end
            cmin2_1(N)=cmin2_1(N)+Penaltyr1+Penaltyr2;

            [cmin2_2(N),A2_1]=min(cmin2_1);
            res_m2_2(N)=res_m2_1(A2_1);
            D_best2_2(N,:)=D_best2_1(A2_1,:);
            T_best2_2(N,:)=T_best2_1(A2_1,:);
            clear cmin2_1;
            clear res_m2_1;
            clear D_best2_1;
            clear T_best2_1;
end

%%stage3
Demand3=Demand*1.03^40;
epanetclose();
epanetloadfile(wdsfile);
% setdata('EN_BASEDEMAND',Demand3);
for N=1:N_poss 
    for n=1:N-1 %start at n=1, since at stage 1 minimum 1 pipe need to be added since without adding any pipe at stage 1, min head requirement is not satisfied
        dd1=D_best2_2(n,:);
        m=(dd1~=0);
        I1=[];
        for q=1:16
            for r=1:16
                if m(q)==0 && q==I(r)
                    I1=[I1 I(r)];
                end
            end
        end
%         [s,s1]=size(I1);
%         p1=factorial(s1)/(factorial(N-n)*factorial(s1-N+n));
%         choices=combntns(I1,N-n);
        choices=randsample(I1,(N-n)*10*(N-n),true);
        choices=reshape(choices,(N-n)*10,(N-n));
        [p1, p2]=size(choices);
%         combos = dec2base((0:14^(N-n)-1)', 14) - '0';
%         combos = sort(combos, 2);
%         combos = unique(combos, 'rows');
        combos=randsample(dia,(N-n)*10*(N-n),true);
        combos=reshape(combos,(N-n)*10,(N-n));
        [s1,s2]=size(combos);
        for i=1:s1
            for j=1:s2
                if combos(i,j)>9
                    combos(i,j)=combos(i,j)-7;
                end
            end
        end
        combos = dia(1 + combos);
        if (N-n)==1
            combos=combos';
        end
        d=zeros(s1,16);
        for i=1:s1
            d(i,:)=dd1;
        end
        for i=1:p1
            cost3=ones(1,s1)*10^40;
            res3_1=zeros(1,s1);
            res3_2=zeros(1,s1);
            res3_3=zeros(1,s1);
%                 d=dd1;
                k=1;
                T(i,:)=T_best2_2(n,:); 
                T1(i,:)=T(i,:);
                for j=1:16
                    if k>(N-n)
                        break;
                    end
                    if j==choices(i,k)
                        for l=1:s1
                            d(l,j)=combos(l,k);
                        end
                        k=k+1;
                        T1(i,j)=40;
                    end
                end
                for j=1:s1
                D=d(j,:);
                m=(T1(i,:)==20);
                for k=1:16
                    if m(k)==1
                        d1(k)=D(k);
                    else
                        d1(k)=0;
                    end
                end
                DD=[Dia d1];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand1);
                ENsolveH();
                h1=getdata('EN_PRESSURE');
        %         pause
                flow1=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h1(l)<Hmin(l)
                        DH1(l)=Hmin(l)-h1(l);
                    else
                        DH1(l)=0;
                    end
                end
                DHmax1=max(DH1);
                Penalty1=DHmax1*10^20;
                tot_h1=h1+ele;
                DD=[Dia dd1];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand2);
                ENsolveH();
                h2=getdata('EN_PRESSURE');
                flow2=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h2(l)<Hmin(l)
                        DH2(l)=Hmin(l)-h2(l);
                    else
                        DH2(l)=0;
                    end
                end
                DHmax2=max(DH2);
                Penalty2=DHmax2*10^20;
                tot_h2=h2+ele;
                DD=[Dia D];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand3);
                ENsolveH();
                h3=getdata('EN_PRESSURE');
                flow3=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h3(l)<Hmin(l)
                        DH3(l)=Hmin(l)-h3(l);
                    else
                        DH3(l)=0;
                    end
                end
                DHmax3=max(DH3);
                Penalty3=DHmax3*10^20;
                tot_h3=h3+ele;
                cost3(j)= TL_lcc_gen_exp(D,T1(i,:),Length)+Penalty1+Penalty2+Penalty3;
                res3_1(j)=Res_TL_exp2(tot_h1,Demand1,flow1);
                if res3_1(j)<0.8
                    Penaltyr1=(0.8-res3_1(j))*10^20;
                else
                    Penaltyr1=0;
                end
                res3_2(j)=Res_TL_exp2(tot_h2,Demand2,flow2);
                if res3_2(j)<0.8
                    Penaltyr2=(0.8-res3_2(j))*10^20;
                else
                    Penaltyr2=0;
                end
                res3_3(j)=Res_TL_exp2(tot_h3,Demand3,flow3);
                if res3_3(j)<0.8
                    Penaltyr3=(0.8-res3_3(j))*10^20;
                else
                    Penaltyr3=0;
                end
                cost3(j)=cost3(j)+Penaltyr1+Penaltyr2+Penaltyr3;
                end
                [cmin3(i),A1]=min(cost3(1,:));
                res_m3(i)=res3_3(A1);
                D_best3(i,:)=d(A1,:);
                T_best3(i,:)=T(i,:);
        end
            [cmin3_1(N),A1]=min(cmin3);
            res_m3_1(N)=res_m3(A1);
            D_best3_1(N,:)=D_best3(A1,:);
            T_best3_1(N,:)=T_best3(A1,:);
            clear cost3;
            clear res3_1;
            clear res3_2;
            clear res3_3;
            clear D;
            clear T;
    end
    D_best3_1(N,:)=D_best2_2(N,:);
    T_best3_1(N,:)=T_best2_2(N,:);
    m=(T_best3_1(N,:)==20);
    for j=1:16
        if m(j)==1
            d1(j)=D_best3_1(j);
        else
            d1(j)=0;
        end
    end
    DD=[Dia d1];
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand1);
    ENsolveH();
    h1=getdata('EN_PRESSURE');
    flow1=getdata('EN_FLOW');
    for l=1:Nnodes
        if h1(l)<Hmin(l)
            DH1(l)=Hmin(l)-h1(l);
        else
            DH1(l)=0;
        end
    end
    DHmax1=max(DH1);
    Penalty1=DHmax1*10^20;
    tot_h1=h1+ele;
    DD=[Dia D_best3_1(N,:)];
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand2);
    ENsolveH();
    h2=getdata('EN_PRESSURE');
    flow2=getdata('EN_FLOW');
    for l=1:Nnodes
        if h2(l)<Hmin(l)
            DH2(l)=Hmin(l)-h2(l);
        else
            DH2(l)=0;
        end
    end
    DHmax2=max(DH2);
    Penalty2=DHmax2*10^20;
    tot_h2=h2+ele;
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand3);
    ENsolveH();
    h3=getdata('EN_PRESSURE');
    flow3=getdata('EN_FLOW');
    for l=1:Nnodes
        if h3(l)<Hmin(l)
            DH3(l)=Hmin(l)-h3(l);
        else
            DH3(l)=0;
        end
    end
    DHmax3=max(DH3);
    Penalty3=DHmax3*10^20;
    tot_h3=h3+ele;
    cmin3_1(N)= TL_lcc_gen_exp(D_best3_1(N,:),T_best3_1(N,:),Length)+Penalty1+Penalty2+Penalty3;
    res1=Res_TL_exp2(tot_h1,Demand1,flow1);
    if res1<0.8
        Penaltyr1=(0.8-res1)*10^20;
    else
        Penaltyr1=0;
    end
    res2=Res_TL_exp2(tot_h2,Demand2,flow2);
    if res2<0.8
        Penaltyr2=(0.8-res2)*10^20;
    else
        Penaltyr2=0;
    end
    res_m3(N)=Res_TL_exp2(tot_h3,Demand3,flow3);
    if res_m3(N)<0.8
        Penaltyr3=(0.8-res_m3(N))*10^20;
    else
        Penaltyr3=0;
    end
    cmin3_1(N)=cmin3_1(N)+Penaltyr1+Penaltyr2+Penaltyr3;
    
    [cmin3_2(N),A3_1]=min(cmin3_1);
    res_m3_2(N)=res_m3(A3_1);
    D_best3_2(N,:)=D_best3_1(A3_1,:);
    T_best3_2(N,:)=T_best3_1(A3_1,:);
    clear cmin3_1;
    clear res_m3_1;
    clear D_best3_1;
    clear T_best3_1;
end

%%stage4
Demand4=Demand*1.03^50;
epanetclose();
epanetloadfile(wdsfile);
% setdata('EN_BASEDEMAND',Demand4);
for N=1:N_poss 
    for n=1:N-1 %start at n=1, since at stage 1 minimum 1 pipe need to be added since without adding any pipe at stage 1, min head requirement is not satisfied
        dd1=D_best3_2(n,:);
        m=(dd1~=0);
        I1=[];
        for q=1:16
            for r=1:16
                if m(q)==0 && q==I(r)
                    I1=[I1 I(r)];
                end
            end
        end
%         [s,s1]=size(I1);
%         p1=factorial(s1)/(factorial(N-n)*factorial(s1-N+n));
%         choices=combntns(I1,N-n);
        choices=randsample(I1,(N-n)*10*(N-n),true);
        choices=reshape(choices,(N-n)*10,(N-n));
        [p1, p2]=size(choices);
        combos = dec2base((0:14^(N-n)-1)', 14) - '0';
        combos = sort(combos, 2);
        combos = unique(combos, 'rows');
        [s1,s2]=size(combos);
        for i=1:s1
            for j=1:s2
                if combos(i,j)>9
                    combos(i,j)=combos(i,j)-7;
                end
            end
        end
        combos = dia(1 + combos);
        if (N-n)==1
            combos=combos';
        end
        d=zeros(s1,16);
        for i=1:s1
            d(i,:)=dd1;
        end
        for i=1:p1
            cost4=ones(1,s1)*10^40;
            res4_1=zeros(1,s1);
            res4_2=zeros(1,s1);
            res4_3=zeros(1,s1);
            res4_4=zeros(1,s1);
%                 d=dd1;
                k=1;
                T(i,:)=T_best3_2(n,:); 
                T1(i,:)=T(i,:);
                for j=1:16
                    if k>(N-n)
                        break;
                    end
                    if j==choices(i,k)
                        for l=1:s1
                            d(l,j)=combos(l,k);
                        end
                        k=k+1;
                        T1(i,j)=50;
                    end
                end
                for j=1:s1
                D=d(j,:);
                m=(T1(i,:)==20);
                for k=1:16
                    if m(k)==1
                        d1(k)=D(k);
                    else
                        d1(k)=0;
                    end
                end
                DD=[Dia d1];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand1);
                ENsolveH();
                h1=getdata('EN_PRESSURE');
        %         pause
                flow1=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h1(l)<Hmin(l)
                        DH1(l)=Hmin(l)-h1(l);
                    else
                        DH1(l)=0;
                    end
                end
                DHmax1=max(DH1);
                Penalty1=DHmax1*10^20;
                tot_h1=h1+ele;
                o=(T1(i,:)==30);
                for k=1:16
                    if m(k)==1 || o(k)==1
                        d1(k)=D(k);
                    else
                        d1(k)=0;
                    end
                end
                DD=[Dia d1];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand2);
                ENsolveH();
                h2=getdata('EN_PRESSURE');
                flow2=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h2(l)<Hmin(l)
                        DH2(l)=Hmin(l)-h2(l);
                    else
                        DH2(l)=0;
                    end
                end
                DHmax2=max(DH2);
                Penalty2=DHmax2*10^20;
                tot_h2=h2+ele;
                p=(T1(i,:)==40);
                for k=1:16
                    if m(k)==1 || o(k)==1 || p(k)==1
                        d1(k)=D(k);
                    else
                        d1(k)=0;
                    end
                end
                DD=[Dia d1];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand3);
                ENsolveH();
                h3=getdata('EN_PRESSURE');
                flow3=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h3(l)<Hmin(l)
                        DH3(l)=Hmin(l)-h3(l);
                    else
                        DH3(l)=0;
                    end
                end
                DHmax3=max(DH3);
                Penalty3=DHmax3*10^20;
                tot_h3=h3+ele;
                DD=[Dia D];
%                 epanetclose();
%                 epanetloadfile(wdsfile);
                setdata('EN_DIAMETER',DD);
                setdata('EN_BASEDEMAND',Demand4);
                ENsolveH();
                h4=getdata('EN_PRESSURE');
                flow4=getdata('EN_FLOW');
                for l=1:Nnodes
                    if h4(l)<Hmin(l)
                        DH4(l)=Hmin(l)-h4(l);
                    else
                        DH4(l)=0;
                    end
                end
                DHmax4=max(DH4);
                Penalty4=DHmax4*10^20;
                tot_h4=h4+ele;
                
                cost4(j)= TL_lcc_gen_exp(D,T1(i,:),Length)+Penalty1+Penalty2+Penalty3+Penalty4;
                res4_1(j)=Res_TL_exp2(tot_h1,Demand1,flow1);
                if res4_1(j)<0.8
                    Penaltyr1=(0.8-res4_1(j))*10^20;
                else
                    Penaltyr1=0;
                end
                res4_2(j)=Res_TL_exp2(tot_h2,Demand2,flow2);
                if res4_2(j)<0.8
                    Penaltyr2=(0.8-res4_2(j))*10^20;
                else
                    Penaltyr2=0;
                end
                res4_3(j)=Res_TL_exp2(tot_h3,Demand3,flow3);
                if res4_3(j)<0.8
                    Penaltyr3=(0.8-res4_3(j))*10^20;
                else
                    Penaltyr3=0;
                end
                res4_4(j)=Res_TL_exp2(tot_h4,Demand4,flow4);
                if res4_4(j)<0.8
                    Penaltyr4=(0.8-res4_4(j))*10^20;
                else
                    Penaltyr4=0;
                end
                cost4(j)=cost4(j)+Penaltyr1+Penaltyr2+Penaltyr3+Penaltyr4;
                end
                [cmin4(i),A1]=min(cost4(1,:));
                res_m4(i)=res4_4(A1);
                D_best4(i,:)=d(A1,:);
                T_best4(i,:)=T(i,:);
        end
            [cmin4_1(N),A1]=min(cmin4);
            res_m4_1(N)=res_m4(A1);
            D_best4_1(N,:)=D_best4(A1,:);
            T_best4_1(N,:)=T_best4(A1,:);
            clear cost4;
            clear res4_1;
            clear res4_2;
            clear res4_3;
            clear res4_4;
            clear D;
            clear T;
    end
    D_best4_1(N,:)=D_best3_2(N,:);
    T_best4_1(N,:)=T_best3_2(N,:);
    m=(T_best4_1(N,:)==20);
    for j=1:16
        if m(j)==1
            d1(j)=D_best4_1(j);
        else
            d1(j)=0;
        end
    end
    DD=[Dia d1];
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand1);
    ENsolveH();
    h1=getdata('EN_PRESSURE');
    flow1=getdata('EN_FLOW');
    for l=1:Nnodes
        if h1(l)<Hmin(l)
            DH1(l)=Hmin(l)-h1(l);
        else
            DH1(l)=0;
        end
    end
    DHmax1=max(DH1);
    Penalty1=DHmax1*10^20;
    tot_h1=h1+ele;
    o=(T_best4_1(N,:)==30);
    for j=1:16
        if m(j)==1 || o(j)==1
            d1(j)=D_best4_1(j);
        else
            d1(j)=0;
        end
    end
    DD=[Dia d1];
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand2);
    ENsolveH();
    h2=getdata('EN_PRESSURE');
    flow2=getdata('EN_FLOW');
    for l=1:Nnodes
        if h2(l)<Hmin(l)
            DH2(l)=Hmin(l)-h2(l);
        else
            DH2(l)=0;
        end
    end
    DHmax2=max(DH2);
    Penalty2=DHmax2*10^20;
    tot_h2=h2+ele;
    p=(T_best4_1(N,:)==40);
    for j=1:16
        if m(j)==1 || o(j)==1 || p(k)==1
            d1(j)=D_best4_1(j);
        else
            d1(j)=0;
        end
    end
    DD=[Dia d1];
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand3);
    ENsolveH();
    h3=getdata('EN_PRESSURE');
    flow3=getdata('EN_FLOW');
    for l=1:Nnodes
        if h3(l)<Hmin(l)
            DH3(l)=Hmin(l)-h3(l);
        else
            DH3(l)=0;
        end
    end
    DHmax3=max(DH3);
    Penalty3=DHmax3*10^20;
    tot_h3=h3+ele;
    DD=[Dia D_best4_1(N,:)];
%     epanetclose();
%     epanetloadfile(wdsfile);
    setdata('EN_DIAMETER',DD);
    setdata('EN_BASEDEMAND',Demand4);
    ENsolveH();
    h4=getdata('EN_PRESSURE');
    flow4=getdata('EN_FLOW');
    for l=1:Nnodes
        if h4(l)<Hmin(l)
            DH4(l)=Hmin(l)-h4(l);
        else
            DH4(l)=0;
        end
    end
    DHmax4=max(DH4);
    Penalty4=DHmax4*10^20;
    tot_h4=h4+ele;
    
    cmin4_1(N)= TL_lcc_gen_exp(D_best4_1(N,:),T_best4_1(N,:),Length)+Penalty1+Penalty2+Penalty3+Penalty4;
    res1=Res_TL_exp2(tot_h1,Demand1,flow1);
    if res1<0.8
        Penaltyr1=(0.8-res1)*10^20;
    else
        Penaltyr1=0;
    end
    res2=Res_TL_exp2(tot_h2,Demand2,flow2);
    if res2<0.8
        Penaltyr2=(0.8-res2)*10^20;
    else
        Penaltyr2=0;
    end
    res3=Res_TL_exp2(tot_h3,Demand3,flow3);
    if res3<0.8
        Penaltyr3=(0.8-res3)*10^20;
    else
        Penaltyr3=0;
    end
    res_m4(N)=Res_TL_exp2(tot_h4,Demand4,flow4);
    if res_m4(N)<0.8
        Penaltyr4=(0.8-res_m4(N))*10^20;
    else
        Penaltyr4=0;
    end
    
    cmin4_1(N)=cmin4_1(N)+Penaltyr1+Penaltyr2+Penaltyr3+Penaltyr4;
    
    [cmin4_2(N),A3_1]=min(cmin4_1);
    res_m4_2(N)=res_m4(A3_1);
    D_best4_2(N,:)=D_best4_1(A3_1,:);
    T_best4_2(N,:)=T_best4_1(A3_1,:);
    clear cmin4_1;
    clear res_m4_1;
    clear D_best4_1;
    clear T_best4_1;
end
[cmin4_3,A4]=min(cmin4_2);
res_m4_3=res_m4_2(A4);
D_best4_3=D_best4_2(A4,:);
T_best4_3=T_best4_2(A4,:);
epanetclose();