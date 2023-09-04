%cleaning up the working environment before starting
clc, clear, close all

%Initiialize Panel Databank
Panels.DB = xlsread('optimizationdatabase','F3:I52');

%Initiialize input parameters
%Assume all panels are square
amb_Temp = 35;
area_Ava = 70;
budget = 20000;
power_Req = 4500;
allowable_Loss_Percent = 15;
allowable_Waste_Percent = 10;
desiredNum = 5;

%Initialize the optimization problem parameters:
Panels.num = 50; %number of panel types
MVO.Num_Generations = 100;  %maximum number of generations to be tested
MVO.Pop_size = 10;  %the population size
MVO.Feasible = 0;
MVO.original_pos = zeros(1,Panels.num);
MVO.WEP_Max=1;
MVO.WEP_Min=0.2;
MVO.Time=1;
Best_universe=zeros(1,dim);
Best_universe_Inflation_rate=inf;
stopping = 0;

MVO.Pop_init = zeros(MVO.Pop_size,Panels.num);  %initialize an empty matrix
MVO.Pop_initB = zeros(MVO.Pop_size,Panels.num);
MVO.Pop_cost = zeros(MVO.Pop_size,1);  %initialize an empty array of zeros

for i = 1:MVO.Pop_size  %for all the population members
    MVO.Feasible = 0;
    while(MVO.Feasible == 0) %loop until sol is feasible
        MVO.Pop_init(i,:) = zeros(1,Panels.num);  %initialize an empty matrix
        MVO.Pop_initB(i,:) = zeros(1,Panels.num);
        MVO.Pop_cost(i,1) = 0;
        cell1 = randi([1,10]);
        MVO.Pop_initB(i,cell1) = 1;
        cell2 = randi([11,20]);
        MVO.Pop_initB(i,cell2) = 1;
        cell3 = randi([21,30]);
        MVO.Pop_initB(i,cell3) = 1;
        cell4 = randi([31,40]);
        MVO.Pop_initB(i,cell4) = 1;
        cell5 = randi([41,50]);
        MVO.Pop_initB(i,cell5) = 1;
        for j = 1:Panels.num  %for all the elements in the array
            val = randi(50);
            MVO.Pop_init(i,j) = val;  %update value of element in solution array
        end

        %Calculate Cost of initial Solution:
        total_Cost = 0;
        total_Power = 0;
        total_Cof= 0;
        total_Num = 0;
        total_Area = 0;
        for j = 1:Panels.num  %for all the elements in the array
            total_Cost = total_Cost + Panels.DB(j,4)*MVO.Pop_init(i,j)*MVO.Pop_initB(i,j); %total cost calculation
            total_Power = total_Power + Panels.DB(j,1)*MVO.Pop_init(i,j)*MVO.Pop_initB(i,j)*0.9; %total energy calculation
            total_Cof = total_Cof + (amb_Temp - 5)*Panels.DB(j,3)*MVO.Pop_init(i,j)*MVO.Pop_initB(i,j); %total temp coeffiecient calculation
            total_Num = total_Num + MVO.Pop_init(i,j)*MVO.Pop_initB(i,j); %total number of panels in current solution calculation
            total_Area = total_Area + MVO.Pop_init(i,j)*Panels.DB(j,2)*MVO.Pop_initB(i,j); %Total area calculation
        end
        area_Waste = (area_Ava-total_Area)/area_Ava; %calculating area wasted percentage
        loss = -(total_Cof/total_Num); %calculating percentage waste due to temp coefficient
        if(total_Cost <= budget && total_Area <= area_Ava && total_Power >= power_Req && loss <= allowable_Loss_Percent && 0 <= area_Waste && area_Waste <= allowable_Waste_Percent) %feasibilty check
            MVO.Feasible = 1; %if feasible set flag to 1
            MVO.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste; %calculating cost of initial sol
            %MVO.Pop_cost(i,1) = total_Cost*1e-4 + (1/total_Power)*1e5 + loss + area_Waste*1e2; %calculating cost of initial sol
        end
    end
end

MVO.Pop_current = MVO.Pop_init;  %save initial population as the current population
MVO.Pop_currentB = MVO.Pop_initB;
clear i j k m cell val total_Area total_Num total_Cof total_Power total_Cost;  %clear used temporary variables

%Plot the cost function figure
Figures.Main_fig = figure;  %create new figure
Best_Cost_Array = zeros(1,MVO.Num_Generations+1);
Best_Cost_Array(1,1) = min(MVO.Pop_cost);  %save the cost of initial solution
figure(Figures.Main_fig);  %access the main figure
subplot(1,2,1);  %in the 1st location in the figure
Figures.Cost = plot(0,Best_Cost_Array(1,1),'m*','LineWidth',1.6); %plot the initial cost
xlim([0,MVO.Num_Generations]); %define limits in X-direction
grid on; %initialize the grid to be on
axis square; %make the axes look like square
xlabel('Number of Generations');
ylabel('Cost function');
title('The Cost function');  %give title to figure
legend('Cost function');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the Multi-Verse Optimizer Loop here   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gen = 1:GW.Num_Generations  %for the maximum number of generations
    if stopping == 200
        break;
    end
    [~, member_index] = sort(GW_new_cost);
    bestSol = GW_new_cost(member_index(1),1);
    Best_Cost_Array(1,gen+1) = bestSol;
    
    WEP = MVO.WEP_Min + gen*((MVO.WEP_Max - MVO.WEP_Min)/GW.Num_Generations);
    TDR=1-((gen)^(1/6)/(GW.Num_Generations)^(1/6));

    Inflation_rates = zeros(1,size(MVO.Pop_current,1));
          
    for j = 1:MVO.Pop_size
        Flag4ub = MVO.Pop_current(gen,:)>ub;
        Flag4lb = MVO.Pop_current(gen,:)<lb;
        MVO.Pop_current(gen,:) = (MVO.Pop_current(gen,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        Inflation_rates(1,i)=fobj((MVO.Pop_current(i,:));
        
        %Elitism
        if Inflation_rates(1,i)<Best_universe_Inflation_rate
            Best_universe_Inflation_rate=Inflation_rates(1,i);
            Best_universe = Universes(i,:);
        end
    end
    [sorted_Inflation_rates,sorted_indexes]=sort(Inflation_rates);
    
    for newindex=1:N
        Sorted_universes(newindex,:)=Universes(sorted_indexes(newindex),:);
    end

    %Normaized inflation rates (NI in Eq. (3.1) in the paper)
    normalized_sorted_Inflation_rates=normr(sorted_Inflation_rates);
    
    MVO.Pop_current(1,:)= Sorted_universes(1,:);
    

end    







