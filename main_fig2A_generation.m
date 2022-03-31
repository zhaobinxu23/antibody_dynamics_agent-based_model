

%% modeling a small group of 1000 people

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
% for i = 1:10000
%     varied_para_1(i) = normrnd(1e-5,0.2e-5) ;
% end
% save('varied_para_1','varied_para_1');
% 
% for i = 1:10000
% varied(i) = min(normrnd(1.8e-8,0.1e-8),1.95e-8);
% end
% save('varied','varied');

varied_para_1 = load('varied_para_1.mat');
varied_para_1 = varied_para_1.varied_para_1;

varied = load('varied.mat');
varied = varied.varied;

Matrix = load('interaction_1000.mat');

    
M = Matrix.interaction_M;

L = Matrix.interaction_M;


delta_t = 1;
virus_complex  = zeros(1000,14000);
virus  =  zeros(1000,14000);
virus(1,1) = 10;
antibody =  100*ones(1000,14000);
environ_antigen =  1e6*ones(1000,14000);
environ_complex =  zeros(1000,14000);
for kk = 1:100
    kk
 
for i = 1:140
    
    for j = 1:1000
        % dd = sum(virus(:,(kk-1)*140+i)'.*M(j,:)/para(10)/(0.9996)^((kk-1)*140+i));
         dd = sum(virus(:,(kk-1)*140+i)'.*M(j,:)/para(10)); %% the transmission capacity of virus is fixed
        if dd <= 0.001
            dd = 0;
        end
        virus_complex(j,(kk-1)*140+i+1) =  max((varied_para_1(j) * virus(j,(kk-1)*140+i) * antibody(j,(kk-1)*140+i) - para(2) * virus_complex(j,(kk-1)*140+i) - para(3) * virus_complex(j,(kk-1)*140+i)) * delta_t + virus_complex(j,(kk-1)*140+i),0);
        antibody(j,(kk-1)*140+i+1) = max((para(2)*virus_complex(j,(kk-1)*140+i) - varied_para_1(j) *antibody(j,(kk-1)*140+i)*virus(j,(kk-1)*140+i) + para(4)*virus_complex(j,(kk-1)*140+i) - para(5)*antibody(j,(kk-1)*140+i) + para(9)*environ_complex(j,(kk-1)*140+i) - varied(j)*antibody(j,(kk-1)*140+i)*environ_antigen(j,(kk-1)*140+i) + para(4)*environ_complex(j,(kk-1)*140+i))*delta_t + antibody(j,(kk-1)*140+i),100);
        virus(j,(kk-1)*140+i+1) = max((-varied_para_1(j) *antibody(j,(kk-1)*140+i)*virus(j,(kk-1)*140+i) + para(2)*virus_complex(j,(kk-1)*140+i) - para(6)*virus(j,(kk-1)*140+i) + (para(7))*virus(j,(kk-1)*140+i) + dd)*delta_t + virus(j,(kk-1)*140+i),0); %% the virus proliferation capacity is also fixed in the ideal model
%         if  virus(j,(kk-1)*140+i+1) > 0.1e5
% %             M = Trace_contact(M,j);
%             M(j,:) = zeros(1,10000);
%         end

        environ_antigen(j,(kk-1)*140+i+1) = 1e6;
        environ_complex(j,(kk-1)*140+i+1) = max((varied(j)*environ_antigen(j,(kk-1)*140+i)*antibody(j,(kk-1)*140+i)-para(9)*environ_complex(j,(kk-1)*140+i)-para(3)*environ_complex(j,(kk-1)*140+i))*delta_t + environ_complex(j,(kk-1)*140+i),0);
        
        
        
    end
end
end




%% plot figure2 A:

severe_case = zeros(1,14000);
for i = 1:14000
    for j = 1:1000
        if virus(j,i) >= 5e4
            severe_case(i) = severe_case(i) + 1;
        end
    end
end

    test_positive = zeros(1,14000);
for i = 1:14000
    for j = 1:1000
        if virus(j,i) >= 0.1e4
            test_positive(i) = test_positive(i) + 1;
        end
    end
end


    asym = zeros(1,14000);
for i = 1:14000
    for j = 1:1000
        if (virus(j,i) >= 0.1e4) && (virus(j,i) <= 0.5e4)
            asym(i) = asym(i) + 1;
        end
    end
end

   mild = zeros(1,14000);
for i = 1:14000
    for j = 1:1000
        if (virus(j,i) > 0.5e4) && (virus(j,i) < 5e4)
            mild(i) = mild(i) + 1;
        end
    end
end
plot(test_positive);
hold on
plot(asym);
hold on
plot(mild);
hold on
plot(severe_case)



%% generate Video 1A


final_matrix = zeros(1000,14000);

for i = 1:14000
    for j = 1:1000
        if virus(j,i) >= 5e4
            final_matrix(j,i) = 1; 
        end
        if (virus(j,i) < 5e4) && (virus(j,i) >=0.5e4)
            final_matrix(j,i) = 0.5; 
        end
        if (virus(j,i) < 0.5e4) && (virus(j,i) >=0.1e4)
            final_matrix(j,i) = 0.25; 
        end

    end
end


dd = load('dd.mat');
dd = dd.dd;


for i = 1:1400
    i
   s = scatter(dd(:,1),dd(:,2),[],1-final_matrix(:,10*i),'*');
   colormap gray
   caxis([0,1]);
%    name = strcat('picture',num2str(i),'.jpg');
%    saveas(gcf,name);

end


