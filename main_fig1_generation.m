
%% modeling 10000 people group, in this group, each person is get infected in the first 1-th time unit. There is no interaction between them. 

clc
clear


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




delta_t = 1; %% represent 0.1day     time unit = 0.1day
virus_complex  = zeros(10000,3000);%% represent 10000 people over 300 days
virus  =  zeros(10000,3000);
virus(:,1) = 10*ones(10000,1); %% each person is infected. 
antibody =  100*ones(10000,3000);
environ_antigen =  1e6*ones(10000,3000);
environ_complex =  zeros(10000,3000);

 
for i = 1:3000
    
    for j = 1:10000
       

        virus_complex(j,i+1) =  max((varied_para_1(j) * virus(j,i) * antibody(j,i) - para(2) * virus_complex(j,i) - para(3) * virus_complex(j,i)) * delta_t + virus_complex(j,i),0);
        antibody(j,i+1) = max((para(2)*virus_complex(j,i) - varied_para_1(j) *antibody(j,i)*virus(j,i) + para(4)*virus_complex(j,i) - para(5)*antibody(j,i) + para(9)*environ_complex(j,i) - varied(j)*antibody(j,i)*environ_antigen(j,i) + para(4)*environ_complex(j,i))*delta_t + antibody(j,i),100);
        virus(j,i+1) = max((-varied_para_1(j) *antibody(j,i)*virus(j,i) + para(2)*virus_complex(j,i) - para(6)*virus(j,i) + (para(7))*virus(j,i))*delta_t + virus(j,i),0);

        environ_antigen(j,i+1) = 1e6;
        environ_complex(j,i+1) = max((varied(j)*environ_antigen(j,i)*antibody(j,i)-para(9)*environ_complex(j,i)-para(3)*environ_complex(j,i))*delta_t + environ_complex(j,i),0);
        
        
        
    end
end

%% plot figure 1A

stdshade(antibody,0.5);

pause



%% plot figure 1B


dd = zeros(10000,1);
for i = 1:10000
    for j = 200:3000
        if (virus(i,j) < virus(i,j-1)) && (virus(i,j) < virus(i,j+1))
            dd(i) = j -  find(virus(i,:)==max(virus(i,:)));
            break
        else 
            dd(i) = 3000 - find(virus(i,:)==max(virus(i,:)));%% dd(i) is the protection duration of individual i. 
        end
    end
end

total = zeros(240,1);
for i = 1:240
    for j = 1:10000
      if dd(j) > 10*i
          total(i) = total(i) + 1;
      end
    end
end


plot(total/100000);