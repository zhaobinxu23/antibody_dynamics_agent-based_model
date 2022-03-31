
clc
clear


%% obtain the morbidity landscape at 30 days when No1 is the origin of this epidemic 
for kkk = 1:1
    kkk
%     clc
%     clear



para(1) = 1e-5; 
para(2) = 1e-14; 
para(3) = 1; 
para(4) = 2; 
para(5) = 0.02; 
para(6) = 0.02; 
para(7) = 0.1; 
% varied(j) = 1.5e-8; 
para(9) = 1e-14;% environment antigen binding kinetics
para(10) = 1e7;% 传播出去的活性病毒总量与体内病毒含量的关系    1e7；


varied_para_1 = load('varied_para_1.mat');
varied_para_1 = varied_para_1.varied_para_1;

varied = load('varied.mat');
varied = varied.varied;

Matrix = load('interaction_1000.mat');

    
M = Matrix.interaction_M;

L = Matrix.interaction_M;


delta_t = 1;
virus_complex  = zeros(1000,400); %% generate the morbidity landscape of those 1000 people over 40 days.
virus  =  zeros(1000,400);
virus(kkk,1) = 10;
antibody =  100*ones(1000,400);
environ_antigen =  1e6*ones(1000,400);
environ_complex =  zeros(1000,400);

 
for i = 1:400
    
    for j = 1:1000
        dd = sum(virus(:,i)'.*M(j,:)/para(10)/(0.9996)^(i));
        if dd <= 0.001
            dd = 0;
        end
        virus_complex(j,i+1) =  max((varied_para_1(j) * virus(j,i) * antibody(j,i) - para(2) * virus_complex(j,i) - para(3) * virus_complex(j,i)) * delta_t + virus_complex(j,i),0);
        antibody(j,i+1) = max((para(2)*virus_complex(j,i) - varied_para_1(j) *antibody(j,i)*virus(j,i) + para(4)*virus_complex(j,i) - para(5)*antibody(j,i) + para(9)*environ_complex(j,i) - varied(j)*antibody(j,i)*environ_antigen(j,i) + para(4)*environ_complex(j,i))*delta_t + antibody(j,i),100);
        virus(j,i+1) = max((-varied_para_1(j) *antibody(j,i)*virus(j,i) + para(2)*virus_complex(j,i) - para(6)*virus(j,i) + (para(7)*0.99995^(i))*virus(j,i) + dd)*delta_t + virus(j,i),0);


        environ_antigen(j,i+1) = 1e6;
        environ_complex(j,i+1) = max((varied(j)*environ_antigen(j,i)*antibody(j,i)-para(9)*environ_complex(j,i)-para(3)*environ_complex(j,i))*delta_t + environ_complex(j,i),0);
        
        
        
    end
end
aa = virus(:,300);

%% transform the virus concentration into symptoms

xx = zeros(1,1000);

for i = 1:1000
    
        if virus(i,300) >= 5e4 %% severe case
            xx(i) = 1; 
        end
        if (virus(j,i) < 5e4) && (virus(j,i) >=0.5e4) %% mild cases
            xx(i) = 0.5; 
        end
        if (virus(j,i) < 0.5e4) && (virus(j,i) >=0.1e4) %% asymptomatic cases
            xx(i) = 0.25; 
        end

    
end

save('current_statue','xx');


end

%% plot figure 6A

dd = load('dd.mat');
dd = dd.dd;  %% this is the geographic information of each individual in this population

s = scatter(dd(:,1),dd(:,2),[],xx,'*');%% the morbidity landscape is displayed, with 1 represents severe case, 0.5 represents mild symptom, 0.25 represents asymptomatic and 0 represents non-infections group. 




