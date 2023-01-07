function [PartclSwormOptimization] = PSO(CostFunction,VarSize,VarMin,VarMax,MaxIt,nPop,w,wdamp,c1,c2,VelMax,VelMin,HamiltonianMemory)
%% Initialization (empty arrays. will contain the particles' values. one array is a sworm)
empty_particle.Position=[]; %position will contain set of rndom values that will get into the cost function 
empty_particle.Cost=[]; %cost will contain the cost function values of the position  
empty_particle.Velocity=[]; %velocity will contain the change that will accure in the prtcl position in the next iteration. depends mainly on the best cost of the whole sworm 
empty_particle.Best.Position=[]; %best position will contain the best position (which gives the lowest cost) of the specific prtcl, acheived through all the iterations so far
empty_particle.Best.Cost=[]; %best position will contain the best cost (lowest) of the specific prtcl, acheived through all the iterations so far
particle=repmat(empty_particle,nPop,1); 
GlobalBest.Cost=inf; %global best will contain the best cost in the whole sworm
for i=1:nPop  %initialization loop. starts the random values that will proceed later in the main loop. each loop forms different prtcl. the whole loop forms the sworm
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize); %formation of npop random positions (one each loop)
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize); %forming a blank vector for velocity
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position); %the value of cost function with the above position
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
end
% Update Global Best
for i=1:nPop
    if particle(i).Best.Cost<GlobalBest.Cost %that if, I think, cant be parallel because it makes a comparison between the whole prtcls in the sworm
        
        GlobalBest=particle(i).Best;
        
    end
    
end
BestCost=zeros(MaxIt,1);
GlobalBest.hamiltonian=zeros(MaxIt,VarSize(2));
%% PSO Main Loop
for it=1:MaxIt %the main loop. each iteration runs on the whole sworm (nPop prtcls)
    GlobalBestPosition=GlobalBest.Position;
    for i=1:nPop %runs on each prtcl
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBestPosition-particle(i).Position);    %velocity. makes a variation of the best position of the specific prtcl and the the best position of the whole sworm and set a change rate for the specific prtcl
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);

        % Evaluation (position value in the cost function)
        particle(i).Cost = CostFunction(particle(i).Position);
        
    end
     % Update Personal Best
     for i=1:nPop %the loop runs on each prtcl, finds personal+global best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost && HamiltonianMemory
                
                GlobalBest.Cost=particle(i).Best.Cost;
                GlobalBest.Position=particle(i).Best.Position;
                GlobalBest.hamiltonian(it,:) = particle(i).Best.Position;
            elseif particle(i).Best.Cost<GlobalBest.Cost

                GlobalBest.Cost=particle(i).Best.Cost;
                GlobalBest.Position=particle(i).Best.Position;
            end
        end
            
     end
     Boolean = sum(GlobalBest.hamiltonian(it,:))==0 && it>2 && sum(GlobalBest.hamiltonian(it-1,:)) ~= 0;
    if Boolean && HamiltonianMemory
        GlobalBest.hamiltonian(it,:) = GlobalBest.hamiltonian(it-1,:);
    end
    BestCost(it)=GlobalBest.Cost; 
    

    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]); %displaying of the iteraion number and best cost

    w=w*wdamp; %updating the inertia weight for the next loop. acts on the velocity
end
%% Results
figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

GlobalBest.memory=BestCost;
PartclSwormOptimization = GlobalBest;

end