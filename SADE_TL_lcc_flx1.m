clear all;
clc;
wdsfile='TL.inp';
addpath('D:\Epanet_codes');
epanetloadfile(wdsfile);
TD=0;
Rmin=0.65;
% tanks = {'14'};
setdata('EN_DURATION',TD);
NP = input('Enter the population size \n');
Imax = input('Enter the maximum number of iterations \n');
Npipes = getdata('EN_LINKCOUNT');
Nnodes = getdata('EN_NODECOUNT')-1;
Length = getdata('EN_LENGTH');
Demand = getdata('EN_BASEDEMAND');
Demand1 =Demand*(1.03)^10;
setdata('EN_BASEDEMAND',Demand1);
ele=getdata('EN_ELEVATION');
n_sim=1000;
dia =[25.4 50.8 76.2 101.6 152.4 203.2 254 304.8 355.6 406.4 457.2 508 558.8 609.6];
Hmin = ones(1,Nnodes)*30;
diaa = randsample(dia,NP*Npipes,true);
disp('Initial Population is \n')
diab=reshape(diaa,NP,Npipes);
% disp('Initial Population is \n')
q=1;
Fm=0.5; CRm=0.5; Fmo(q)=Fm; CRmo(q)=CRm;
netcost = zeros(Imax,NP);
initialcost=zeros(1,NP);
head=zeros(1,Nnodes);
index=zeros(1,NP);
% Rel =zeros(Imax,NP);
% n_rel=zeros(NP,Nnodes);
% index1=zeros(1,NP);
% initialrel=zeros(1,NP);
Best=ones(1,Imax);
DH=zeros(1,Nnodes);
for i=1:NP
    diac =diab(i,:);
    setdata('EN_DIAMETER',diac);
    ENsolveH();
    head =getdata('EN_PRESSURE');
    flow=getdata('EN_FLOW');
    tot_h=head+ele;
    for j=1:Nnodes
        if(head(j)<Hmin(j))
            DH(j)=Hmin(j)-head(j);
        else
            DH(j)=0;
        end
    end
    DHmax=max(DH);
    Penalty=DHmax*10^9;
    initialcost(i)=TL_lcc(diac,Length)+Penalty;
    initialrel(i)=Res_TL(tot_h,Demand1,flow);
    if initialrel(i)<Rmin
        Penaltyr=(Rmin-initialrel(i))*10^10;
    else
        Penaltyr=0;
    end
    initialcost(i)=initialcost(i)+Penaltyr;
end
% disp('Initial reliability values are \n')
% initialrel
% disp('Initial cost is \n')
initialcost;
% wdsfile='mincut.inp';
% epanetloadfile(wdsfile);
Best1=1;
for i=1:NP
    if(initialcost(i)<initialcost(Best1))
        Best1=i;
    end
end
F =abs(randn(1,NP)*0.3+Fm);
CR = abs(randn(1,NP)*0.1+CRm);
% strategy =rand(NP,Npipes);
comp = rand(NP,Npipes);
v = zeros(NP,Npipes);
u = zeros(1,Npipes);
ud = zeros(NP,Npipes);
sf=zeros(Imax,NP);
% sc=zeros(Imax,Npipes);
sflag=0; it=1;
tic
while (it<=Imax && sflag <=200)
    netcost1=zeros(1,NP);
%     Rel1=zeros(1,NP);
    index = zeros(1,NP);
%     index1=zeros(1,NP);
    for i=1:NP
        r = randi(NP,1,3);
        if(i==NP)
            for j=1:3
                if(r(j)==i)
                    r(j)=i-1;
                end
            end
        else
            for j=1:3
                if(r(j)==i)
                    r(j)=i+1;
                end
            end
        end
        for j=1:Npipes
            v(i,j)= diab(r(1),j)+ F(1,i)*(diab(r(2),j)-diab(r(3),j));
        end
        for j=1:Npipes
            if(comp(i,j)<=CR(1,i))
                u(j)=v(i,j);
%                 sc(it,j)=1;
            else
                u(j) = diab(i,j);
            end
        end
%         wdsfile='mincut.inp';
%         epanetloadfile(wdsfile);
        ud(i,:)=Discrete_TL(u);
        setdata('EN_DIAMETER',ud(i,:));
        ENsolveH();
        head = getdata('EN_PRESSURE');
        flow=getdata('EN_FLOW');
        tot_h=head+ele;
        for j=1:Nnodes
            if(head(j)<Hmin(j))
                DH(j)=Hmin(j)-head(j);
            else
                DH(j)=0;
            end
        end
        DHmax=max(DH);
        Penalty=DHmax*10^9;
       netcost1(i)=TL_lcc(ud(i,:),Length)+Penalty;
        Rel1(i)=Res_TL(tot_h,Demand1,flow);
        if Rel1(i)<Rmin
            Penaltyr=(Rmin-Rel1(i))*10^10;
        else
            Penaltyr=0;
        end
        netcost1(i)=netcost1(i)+Penaltyr;
        netcost1;
%         Rel1;
      if(it==1)
        if((netcost1(i)<initialcost(i)))% && Rel1(i)>=0.90)
            Rel(it,i)=Rel1(i);
            netcost(it,i)=netcost1(i);
            diab(i,:)=ud(i,:);
            sf(it,i)=1;
        else
            Rel(it,i)=initialrel(i);
            netcost(it,i)=initialcost(i);
        end
   else
    if((netcost1(i)<netcost(it-1,i)))% && Rel1(i)>=0.90)
        Rel(it,i)=Rel1(i);
        netcost(it,i)=netcost1(i);
        diab(i,:)=ud(i,:);
        sf(it,i)=1;
    else
        Rel(it,i)=Rel(it-1,i);
        netcost(it,i)=netcost(it-1,i);
    end
  end
    end  
 for i=1:NP
    if(netcost(it,i)<netcost(it,Best(it)))
        Best(it)=i;
    end
 end
  if(mod(it,10)==0)
      for x=(it-9):it
      a=0;b=0;cf=0;
      for i=1:NP
         if(sf(it,i)==1)
           a=a+F(i);
           b=b+CR(i);
           cf=cf+1;
         end
      end
      end
      if(cf~=0)
        Fm=a/cf;
        CRm=b/cf;
      end
      q=q+1;
      Fmo(q)=Fm; CRmo(q)=CRm;
     F =abs(randn(1,NP)*0.3+Fm);
    CR = abs(randn(1,NP)*0.1+CRm);
    
  end
  if it>1
      if abs(netcost(it,Best(it))-netcost(it-1,Best(it-1)))<0.0001
          sflag=sflag+1;
      else
          sflag=0;
      end
  end
  it=it+1;
end
% Rel
netcost
diab
Fm
CRm
toc
% wdsfile='mincut.inp';
% epanetloadfile(wdsfile);
epanetclose();
bestcost=0;
for i=1:it-1
    bestcost(i)=netcost(i,Best(i));
end
bestcost(it-1)
diab(Best(it-1),:)
Rel(it-1,Best(it-1))