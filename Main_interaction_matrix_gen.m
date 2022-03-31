clc
clear

%% generate correct interaction_matrix 

clc
clear

%% 4 towns , ff: town 1，gg:town 2，ee: town 3，hh: town4，each town have 250 people
ff = 100/(10^0.5)*rand(250,2);
gg = 100/(10^0.5)*rand(250,2)+41.6228;
ee = 100/(10^0.5)*rand(250,2)+41.6228*2;
hh = 100/(10^0.5)*rand(250,2)+41.6228*3;
dd = [ff;gg;ee;hh];



distance = zeros(1000,1000);
interaction_M = zeros(1000,1000);
for i = 1: 1000
    for j  = 1:1000
        distance(i,j) = (dd(i,1) -dd(j,1))^2 + (dd(i,2) -dd(j,2))^2;
    end
end
for i = 1: 1000
    for j  = 1:1000
        if (distance(i,j) == 0) || (distance(i,j) > 25) %% no contact if the distance is bigger than 5
            interaction_M(i,j) = 0;
        
        else


              interaction_M(i,j) = min(0.8,10/distance(i,j)^2);
        end

    
    end
end




%% add some interaction linkage among those 4 towns
interaction_M(250,251) = 0.8;
interaction_M(251,250) = 0.8;
interaction_M(500,501) = 0.8;
interaction_M(501,500) = 0.8;
interaction_M(750,751) = 0.8;
interaction_M(751,750) = 0.8;
interaction_M(1000,1) = 0.8;


save('interaction_1000','interaction_M');

save('dd','dd');
