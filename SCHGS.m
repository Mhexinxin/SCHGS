function [bestPositions,Convergence_curve]=HGS(N,Max_iter,lb,ub,dim,fobj)
disp('HGS_FEs is now tackling your problem')
p2=0.1;
%tic
% initialize position
bestPositions=zeros(1,dim); %�?优�?�（只有�?个最优�?�，只有�?行）
tempPosition=zeros(N,dim);  %存放每个个体的�??

Destination_fitness=inf;%change this to -inf for maximization problems
Worstest_fitness=-inf;
AllFitness = inf*ones(N,1);%record the fitness of all positions N�?1�?
VC1 = ones(N,1);%record the variation control of all positions

weight3 = ones(N,dim);%hungry weight of each position
weight4 = ones(N,dim);%hungry weight of each position

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=[];

it=1; %Number of iterations
FEs=0;   %Number of current evaluations

hungry = zeros(1,size(X,1));%record the hungry of all positions
count=0;

% Main loop
while  FEs < Max_iter
    VC2 = 0.03; %The variable of variation control 
    
    sumHungry = 0;%record the sum of each hungry 
    
    %sort the fitness
    for i=1:size(X,1)
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
        FEs = FEs+1;

    end
    %循环完之后，可以得到�?好的适应度�?�和�?坏的适应度�??
    
    [AllFitnessSorted,IndexSorted] = sort(AllFitness);
    bestFitness = AllFitnessSorted(1);%按列遍历，第�?个就是最好的
    worstFitness = AllFitnessSorted(size(X,1));%�?后一个是�?坏的
    
    %利用�?优的适应度�?�找出最好的位置（每次迭代先找出�?个，后面再更新）
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(IndexSorted(1),:);
        Destination_fitness = bestFitness;
        count=0;
    end
    
    if worstFitness > Worstest_fitness
        Worstest_fitness = worstFitness;
    end
    
    %计算每个位置的饥饿�?�和总饥饿�??
    for i = 1:size(X,1)
         %calculate the variation control of all positions
         VC1(i) = sech(abs(AllFitness(i)-Destination_fitness));    
         %calculate the hungry of each position
        if Destination_fitness == AllFitness(i) %如果满足�?优的适应度�??
            hungry(1,i) = 0;
            count = count+1;
            tempPosition(count,:)=X(i,:); %将局部最优的位置存起�?
        else
            temprand = rand();
            c = (AllFitness(i)-Destination_fitness)/(Worstest_fitness-Destination_fitness)*temprand*2*(ub-lb);
            if c<100
                b=100*(1+temprand);
            else
                b=c;
            end   
            hungry(1,i) = hungry(1,i)+ max(b); 
            sumHungry = sumHungry + hungry(1,i);
        end
    end 
    
    %calculate the hungry weight of each position 计算饥饿权重
    for i=1:size(X,1)
        for j=2:size(X,2)
                weight3(i,j) = (1-exp(-abs(hungry(1,i)-sumHungry)))*rand()*2;
                if rand()<VC2
                    weight4(i,j) = hungry(1,i)*size(X,1)/sumHungry*rand();
                else
                    weight4(i,j) = 1;
                end
        end
        
    end
    
    
    % Update the Position of search agents
    shrink=2*(1-FEs/Max_iter); % a decreases linearly fron 2 to 0
    Xnew=X;
    
    for i=1:size(X,1)
        if rand<VC2  %探索阶段
            Xnew(i,:) = X(i,j)*(1+randn(1));
        else  %�?发阶�?
            A = randi([1,count]);
            for j=1:size(X,2)
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    Xnew(i,j) = weight4(i,j)*tempPosition(A,j)+vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));
                else
                    Xnew(i,j) = weight4(i,j)*tempPosition(A,j)-vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));
                end
            end
        end
    end
    CW=size(Xnew)
   for i=1:N
        %Boundary absorption
        Flag4ub=Xnew(i,:)>ub;
        Flag4lb=Xnew(i,:)<lb;
        Xnew(i,:)=(Xnew(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newRime_rates(i)=fobj(Xnew(i,:));
        FEs=FEs+1;
        %Positive greedy selection mechanism
        if newRime_rates(i)<AllFitness(i)
            AllFitness(i) = newRime_rates(i);
            X(i,:) = Xnew(i,:);
            if newRime_rates(i)< Destination_fitness
               Destination_fitness=AllFitness(i);
               bestPositions=Xnew(i,:);
            end
        end
   end
    
   %%  CrissCross
    [newX,fitness,FEs] = CC_p(X,AllFitness,dim,lb,ub,fobj,FEs,p2);
    for i=1:N
      if fitness(i)<AllFitness(i)
         X(i,:)=newX(i,:);
         AllFitness(i)=fitness(i);
      end 
    end
   %%
    
    
    
    Convergence_curve(it)=Destination_fitness;
    it=it+1;
end
%toc
end




