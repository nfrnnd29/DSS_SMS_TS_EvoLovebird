clc
clear

%% Dataset
Job_processing_time = 5; %3, 5,2
number_instance = 14;%10, 5,0
number_jobs = 10;
%instance_str = ( [num2str(Job_processing_time) '_' num2str(number_instance) '_' num2str(number_jobs) '.txt']);
instance_name = ([num2str(number_instance) '-' num2str(number_jobs) 'x' num2str(Job_processing_time)]);
instance_extension = ([num2str(number_instance) '_' num2str(number_jobs) 'x' num2str(Job_processing_time) '.txt']);
instance = importdata(instance_extension);

%%
%instance = importdata('3_8.txt');%0_10x3_instance.txt, 0_10x5_instance.txt, 12_10x3_instance.txt
MaxGen = 1000; %500,1000,8000
agent = 10000; %100,500,4000
n_job = instance(1,1);
instance = instance';
instance = instance(2:end);
a_all = {};
instances  = instance;

for c = 1:3
    a1 = instances(1:number_jobs); %(1:(length(instance)/instance(1,1))) ; instances(1:(length(instance)/instance(1,1)))
    instances(1:number_jobs)=[];%instances(1:(length(instance)/instance(1,1))) = [];
    a_all(end+1) = {a1};
end

processing_time = cell2mat(a_all(1,1));
p_time = processing_time(1,1);
release_date_all = cell2mat(a_all(1,2))';
weight_job_all = cell2mat(a_all(1,3))';

P_idx = 1:n_job;

%% Index Job Table
n_square_mat = n_job*p_time;
[AA,full_weight_mat,indexed_fix_before_summation,minjobparts,maxjobparts] = assign_mat_p2(n_job,p_time,release_date_all,weight_job_all);

%% Search Minimum Release Date
assign_table = indexed_fix_before_summation;
odd_index = 1:p_time:size(indexed_fix_before_summation,1);
first_job_part = indexed_fix_before_summation(odd_index,:);
placed_schedule = {};
as_placed_mat = {};
compare_weight = {};
compare_weight_idx = {};
weight_best = {};

while ~isnan(first_job_part.SchedulePlacement(:)==p_time*n_job)
    [value, index] = min(first_job_part.SchedulePlacement(:));
    if isnan(value)
        break
    end
    [value_same, index_same] = find(first_job_part.SchedulePlacement == value);
    if length(index_same) ~= 1 %jika ukurannya lebih dari satu
        for g = 1:length(index_same) %menemukan minimum lainnya
            value_same(g);
            value_s1 = first_job_part.JobWeight(value_same(g));
            index_s1 = first_job_part.MatrixIndex(value_same(g));
            compare_weight{end+1} = value_s1; %dimasukkan variabel
            compare_weight_idx{end+1} = index_s1;
        end
        compare_fix = max(cell2mat(compare_weight));
        value_idx_f = find(cell2mat(compare_weight) == compare_fix);
        value_fix = value_same(value_idx_f);

        if length(value_idx_f) ~= 1
            placed_schedule(end+1) = {first_job_part.JobIndex(value_fix(1,1))};
            as_placed_mat(end+1) = {first_job_part.MatrixIndex(value_fix(1,1))+1};
            first_job_part.SchedulePlacement(value_fix(1,1)) = NaN;

        elseif length(value_idx_f) == 1
            placed_schedule(end+1) = {first_job_part.JobIndex(value_fix)};
            as_placed_mat(end+1) = {first_job_part.MatrixIndex(value_fix)+1};
            first_job_part.SchedulePlacement(value_fix) = NaN;
        end
        compare_weight = {};
        compare_weight_idx = {};
        value_fix = 0;
        compare_fix = 0;
        value_idx_f = 0;
    else
        placed_schedule(end+1) = {first_job_part.JobIndex(index)};
        as_placed_mat(end+1) = {first_job_part.MatrixIndex(index)+1};
        first_job_part.SchedulePlacement(index) = NaN;
    end
    [value_one, index_one] = min(first_job_part.SchedulePlacement(:));
    first_job_part.SchedulePlacement(index_one) = first_job_part.SchedulePlacement(index_one) + 1;
end


as_placed_mat;
placed_schedule;

as_placed = cell2mat(placed_schedule(1));

%% Create Initial Solution
initial = zeros(agent,n_square_mat);
b = repelem(1:n_job,p_time);
[uv_initial,uv_initial_idx,~] = unique(b,'first','legacy');
%[~,find_first_idx] = find(uv_initial == as_placed);
b(uv_initial_idx(as_placed)) = [];
%v = b(2:end);
for p = 1:agent
    initial(p,:) = [as_placed b(randperm(length(b)))];
end
initial;


%% Calculate Objective Function
%"0,0,2,2,4,4,1,1,9,9,8,8,3,3,5,5,7,7,6,6"
%1 3 5 2 10 9 4 6 8 7 ==> 1484 ==> BEST FITNESS !!!!!


evenI = 2:2:size(indexed_fix_before_summation,1);
oddI = 1:2:size(indexed_fix_before_summation,1);
lastI = 1:p_time:size(indexed_fix_before_summation,1);

% reorder matrix to reconstruct original matrix
%second_job_part = indexed_fix_before_summation(evenI,:);
second_job_part = indexed_fix_before_summation(lastI,:);

second_job_part_s = second_job_part(:,[1 2 3 5 6]); %1 2 5]
ShuffleIndex = (second_job_part.SchedulePlacement) * 0;
ConstraintShuffle = (second_job_part.SchedulePlacement) * 0;
second_job_part_s = [second_job_part_s array2table(ShuffleIndex) array2table(ConstraintShuffle)];


second_job_part_s;
check_initial = initial;
second_job_part_shuffle = second_job_part_s;

c = 0;
pinalty = 9999999;
cost = NaN(agent,1);
for i = 1:agent
    n_initial_solution = i;
    intial_job_check = check_initial(i,:);
    %intial_job_check = [1 1 3 3 5 5 2 2 10 10 9 9 4 4 6 6 8 8 7 7]
    [unique_value,first_jp_idx,~] = unique(intial_job_check,'first','legacy');
    [~,second_jp_idx,~] = unique(intial_job_check,'last','legacy');
    %first job part check
    for j = 1:n_job
        job_part_check = unique_value(j);
        time_release = second_job_part_shuffle.JobReleaseTime(job_part_check);
        if first_jp_idx(j) >= time_release
            %second job check
            %job_constraint = second_job_part_shuffle.ShuffleIndex(job_part_check)
            job_weight = second_job_part_shuffle.JobWeight(job_part_check);
            if second_jp_idx(j) >= time_release+(p_time-1) %1 %2*job_constraint
                c = c + second_jp_idx(j)*job_weight - time_release;
                cost(i) = c;
            else
                cost(i) = pinalty;
                break
            end
        else
            cost(i) = pinalty;
            break
        end
    end
    c = 0;
end

pbest = check_initial;
f_pbest = cost;

[f_gbest,ind_f] = min(f_pbest);
gbest = pbest(ind_f,:);
BestFitIter(1) = f_gbest;

% %% Sorted schedule
% v = {};
% for i = 1:length(placed_schedule)
%     number = placed_schedule(i);
%     for x = 1:p_time
%         v{end+1} = cell2mat(number);
%     end
% end
% 
% v =cell2mat(v);
% 
% %% Sorted schedule Cost
% 
% c = 0;
% v_cost = NaN;
% intial_job_check = v;
% %intial_job_check = [1 1 3 3 5 5 2 2 10 10 9 9 4 4 6 6 8 8 7 7]
% [unique_value,first_jp_idx,~] = unique(intial_job_check,'first','legacy');
% [~,second_jp_idx,jp] = unique(intial_job_check,'last','legacy');
% %first job part check
% for j = 1:n_job
%     job_part_check = unique_value(j);
%     time_release = second_job_part_shuffle.JobReleaseTime(job_part_check);
%     if first_jp_idx(j) >= time_release
%         %second job check
%         %job_constraint = second_job_part_shuffle.ShuffleIndex(job_part_check)
%         job_weight = second_job_part_shuffle.JobWeight(job_part_check);
%         if second_jp_idx(j) >= time_release+(p_time-1) %1 %2*job_constraint
%             c = c + second_jp_idx(j)*job_weight;
%             v_cost = c;
%         else
%             pinalty = 9999999;
%             v_cost = pinalty;
%             break
%         end
%     else
%         pinalty = 9999999;
%         v_cost = pinalty;
%         break
%     end
% end


%% Problem Definition

%model = CreateModel();      % Create TSP Model

ActionList=CreatePermActionList(p_time*n_job);%CreatePermActionList(2*n_job);    % Action List

nAction=numel(ActionList);              % Number of Actions


%% Tabu Search Parameters

MaxIt=1000; %MaxGen;                      % Maximum Number of Iterations

TL=round(0.5*nAction); %round(0.5*nAction);      % Tabu Length, 0.25 bisa cepet ketemu

%% Initialization

% Create Empty Individual Structure
empty_individual.Position=[];
empty_individual.Cost=[];

% Array to Hold Best Costs
BestCost=zeros(MaxIt,1);
%BestSchedule =zeros(MaxIt,2*n_job);
BestSchedule = {}';

% Initialize Action Tabu Counters
TC=zeros(nAction,1);


%% Initialization Tabu Search

% Create Initial Solution
sol=empty_individual;
sol.Position=gbest;%v;
sol.Cost=f_gbest;%v_cost;

% Initialize Best Solution Ever Found
BestSol=sol;

%% Tabu Search Main Loop
tic;
tStart = cputime;
for it=1:MaxIt %MaxGen

    bestnewsol.Cost=inf;

    % Apply Actions
    for i=1:nAction
        if TC(i)==0
            newsol.Position=DoAction(sol.Position,ActionList{i});
            [unique_value,first_jp_idx,~] = unique(newsol.Position,'first','legacy');
            [~,second_jp_idx,jp] = unique(newsol.Position,'last','legacy');
            %first job part check
            for j = 1:n_job
                job_part_check = unique_value(j);
                time_release = second_job_part_shuffle.JobReleaseTime(job_part_check);
                if first_jp_idx(j) >= time_release
                    %second job check
                    %job_constraint = second_job_part_shuffle.ShuffleIndex(job_part_check)
                    job_weight = second_job_part_shuffle.JobWeight(job_part_check);
                    if second_jp_idx(j) >= time_release+(p_time-1) %1 %2*job_constraint
                        c = c + second_jp_idx(j)*job_weight - time_release;
                        costs = c;
                    else
                        pinalty = 9999999;
                        costs = pinalty;
                        break
                    end
                else
                    pinalty = 9999999;
                    costs = pinalty;
                    break
                end
            end
            c = 0;
            newsol.Cost=costs;
            newsol.ActionIndex=i;

            if newsol.Cost<=bestnewsol.Cost
                bestnewsol=newsol;
            end
        end
    end

    % Update Current Solution
    sol=bestnewsol;

    % Update Tabu List
    for i=1:nAction
        if i==bestnewsol.ActionIndex
            TC(i)=TL;               % Add To Tabu List
        else
            TC(i)=max(TC(i)-1,0);   % Reduce Tabu Counter, max
        end
    end

    % Update Best Solution Ever Found
    if sol.Cost<=BestSol.Cost
        BestSol=sol;
    end

    % Save Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BestSchedule(it) = mat2cell(BestSol.Position,1);
    %abc = cell2mat(BestSol.Position,1);
    abc = cell2mat(BestSchedule(it));

    %BestSchedule = BestSchedule';
    %f_gbest = min(BestCost);
    %[f_pbestt,ind_f] = min(BestCost);
    %gbest = cell2mat(BestSchedule); %pbest(ind_f,:);


    BestFitIter(it) = BestCost(it);

    %best_schedule = gbest;

    % Show Iteration Information
    disp(['Iteration of TS ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ' Best Schedule =' num2str(abc)]);
    %disp(['Iteration of TS ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);% ' Best Schedule =' num2str(abc)]);
end
%disp(['Iteration of TS ' num2str(Gen) ': Best fitness =' num2str(BestFitIter(Gen)) ' Best Schedule =' num2str(best_schedule)]);

toc

%% Plot
plot(1:MaxGen,BestFitIter);
xlabel('Iteration');
ylabel('Best fitness');
title('1-10x2 Dataset Instance with Tabu Search - Lovebird Algorithm.');