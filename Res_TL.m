function resiliency = Res_TL(head,dem,flow)
Hmin=[180 190 185 180 195 190];
t1=0;
t2=0;
ND=6; 
for i=1:ND
    t1=t1+(dem(i)*(head(i)));
    t2=t2+(dem(i)*Hmin(i));
end
t3=flow(1)*210;
resiliency=(t1-t2)/(t3-t2);
