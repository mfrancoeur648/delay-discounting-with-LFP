%DD Behavior Plotting
%% Low Dose 
R128_5sn2 = [2.43 1.79]
R128_5sn2 = [70.9 64.2]
R128_20sn2 = [78.2 1]
R161_5sn2 = [76.9 96.9]
R161_20sn2 = [81.48 4/13]
R161_20sn2 = [81.48 100*4/13]
R163_5sn2 = [100*96/106 100]
R143_5sn2 = [100*8/20 100*18/24]
R143_20sn2 = [100*2/30 100*5/11]
R165_5sn2 = [100*116/117 100*79/80]
R136_5sn2 = [100*3/4 100*2/5]
R164_5sn2 = [100*165/167 100*134/136]
R164_20sn2 = [100*50/52]
R165_20sn2 = [100*24/27]
R163_20sn2 = [100*21/35 100*22/24]
mat_5sn2 = [R128_5sn2' R136_5sn2' R143_5sn2' R161_5sn2' R163_5sn2' R164_5sn2' R165_5sn2']
mean_5sn2 = mean(mat_5sn2)
mat_5sn2(1)
mat_5sn2(1,:)
salavg_5sn2 = mean(mat_5sn2(1,:))
ketavg_5sn2 = mean(mat_5sn2(2,:))
R164_20sn2 = [R164_20sn2 1]
R165_20sn2 = [R165_20sn2 1]
mat_5sn2 = [R128_5sn2' R136_5sn2' R143_5sn2' R161_5sn2' R163_5sn2' R164_5sn2' R165_5sn2']
mat_20sn2 = [R128_20sn2' R143_20sn2' R161_20sn2' R163_20sn2' ]
salavg_20sn2 = mean(mat_20sn2(1,:))
ketavg_20sn2 = mean(mat_20sn2(2,:))
saline = [salavg_5sn2 salavg_20sn2]
ketamine = [ketavg_5sn2 ketavg_20sn2]

figure;
bar([saline,ketamine]);



%% Medium Dose 
R136_M5 = 100*[59/61 40/58]; 
R128_M5 = 100*[NaN 89/170];
R161_M5 = 100*[NaN 1];
R143_M5 = 100*[31/42 71/75]; 
R164_M5 = 100*[141/144 107/112]; 
R165_M5 = 100*[86/87 107/111]; 
R163_M5 = 100*[63/64 112/114]; 

R136_M20 = 100*[20/28 16/19]; 
R128_M20 = 100*[48/106 24/79];
R161_M20 = 100*[43/43 NaN];
R143_M20 = 100*[11/13 NaN]; 
R164_M20 = 100*[70/71 30/32]; 
R165_M20 = 100*[43/43 4/4]; 
R163_M20 = 100*[71/75 13/13]; 

med_5sn2 = [R128_M5' R136_M5' R143_M5' R161_M5' R163_M5' R164_M5' R165_M5']
med_20sn2 = [R128_M20' R136_M20' R143_M20' R163_M20' R164_M20' R165_M20']

ketmed_5sn2 = nanmean(med_5sn2(1,:))
salmed_5sn2 = nanmean(med_5sn2(2,:))

ketmed_20sn2 = nanmean(med_20sn2(1,:))
salmed_20sn2 = nanmean(med_20sn2(2,:))

sal_med = [salmed_5sn2 salmed_20sn2]
ket_med = [ketmed_5sn2 ketmed_20sn2]
figure;
bar([sal_med,ket_med]);

%% No Drug (post ketamine) 
R136_nd20 = 100*[1/22 2/36 5/19 0/11]; 
R128_nd20 = 100*[24/26 NaN 44/85 2/4];
R161_nd20 = 100*[51/67 37/37 53/53 19/20];
R143_nd20 = 100*[3/9 71/91 99/102 97/102]; 
R164_nd20 = 100*[73/150 74/75 66/72 18/21]; 
R163_nd20 = 100*[58/77 71/91 57/89 40/50]; 
R165_nd20 = 100*[NaN 46/48 NaN 42/42]; 

  
R136_nd5 = 100*[3/41 17/46 NaN]; 
R128_nd5 = 100*[NaN 192/198 75/79];
R161_nd5 = 100*[17/17 24/25 NaN];
R143_nd5 = 100*[72/83 71/90 176/182]; 
R164_nd5 = 100*[111/129 136/145 78/82]; 
R163_nd5 = 100*[72/99 78/102 16/25]; 
R165_nd5 = 100*[117/117 70/70 NaN];   



ND_5 = [R128_nd5' R136_nd5' R143_nd5' R161_nd5' R163_nd5' R164_nd5' R165_nd5']
ND_20 = [R128_nd20' R136_nd20' R143_nd20' R161_nd20' R163_nd20' R164_nd20' R165_nd20']

labels = ["R128" "R136" "R143" "R161" "R163" "R164" "R165"]; 
S = size(ND_5)
 
for i = 1:S(2) 
    nd5(i) = nanmean(ND_5(:,i));
    nd20(i) = nanmean(ND_20(:,i));
    figure;
   % bar([nd5,nd20]);
    title(labels(i));  
end

nd5_bar = nanmean(nd5);
nd20_bar = nanmean(nd20);
figure;
bar([nd5_bar,nd20_bar])


%% Pre Drug Behvioral Data

R136_pd5 = 100*[77/90 73/73 NaN NaN]; 
R128_pd5 = 100*[NaN NaN NaN NaN];
R161_pd5 = 100*[55/60 48/59 NaN NaN];
R143_pd5 = 100*[2/74 NaN NaN NaN]; 
R164_pd5 = 100*[157/162 NaN NaN NaN]; 
R163_pd5 = 100*[NaN NaN NaN NaN]; 
R165_pd5 = 100*[2/5 NaN NaN NaN]; 


R136_pd10 = 100*[53/53 28/34 18/24 30/32]; 
R128_pd10 = 100*[99/112 NaN NaN NaN];
R161_pd10 = 100*[NaN NaN NaN NaN];
R143_pd10 = 100*[16/42 16/48 1/26 18/34]; 
R164_pd10 = 100*[NaN 61/64 74/77 95/99]; 
R163_pd10 = 100*[NaN NaN NaN 33/46]; 
R165_pd10 = 100*[NaN NaN NaN 33/46];   


R136_pd20 = 100*[NaN NaN NaN NaN]; 
R128_pd20 = 100*[49/74 45/52 NaN NaN];
R161_pd20 = 100*[19/30 9/22 NaN NaN];
R143_pd20 = 100*[72/83 71/90 176/182 NaN]; 
R164_pd20 = 100*[NaN NaN NaN NaN]; 
R163_pd20 = 100*[35/41 54/56 NaN NaN]; 
R165_pd20 = 100*[6/23 NaN NaN NaN];   
  

R136_pd40 = 100*[3/4 NaN NaN NaN]; 
R128_pd40 = 100*[26/45 27/31 NaN NaN];
R161_pd40 = 100*[6/6 6/14 14/26 NaN];
R143_pd40 = 100*[21/27 NaN NaN NaN]; 
R164_pd40 = 100*[61/70 NaN NaN NaN]; 
R163_pd40 = 100*[27/36 NaN NaN NaN]; 
R165_pd40 = 100*[19/28 4/19 NaN NaN];   

pd_5 = [R128_pd5' R136_pd5' R143_pd5' R161_pd5' R163_pd5' R164_pd5' R165_pd5']
pd_20 = [R128_pd20' R136_pd20' R143_pd20' R161_pd20' R163_pd20' R164_pd20' R165_pd20']
pd_10 = [R128_pd10' R136_pd10' R143_pd10' R161_pd10' R163_pd10' R164_pd10' R165_pd10']
pd_40 = [R128_pd40' R136_pd40' R143_pd40' R161_pd40' R163_pd40' R164_pd40' R165_pd40']

labels = ["R128" "R136" "R143" "R161" "R163" "R164" "R165"]; 
S2 = size(pd_5)
 
for i = 1:S2(2) 
    pd5 = nanmean(pd_5(:,i));
    pd20 = nanmean(pd_20(:,i));
    pd10 = nanmean(pd_10(:,i));
    pd40 = nanmean(pd_40(:,i));
    figure;
    bar([pd5,pd10,pd20,pd40]);
    title(labels(i));  
end



