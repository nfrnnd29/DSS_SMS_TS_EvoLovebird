function [assign_mat,full_weight_mat,indexed_fix_before_summation, minjobparts,maxjobparts]  = assign_mat_p2(n_job, p_time,release_date_all,weight_job_all)
%assign_mat = zeros(n_square_mat,n_square_mat)

%assign_mat_p1(instance, n_square_mat,P_idx,job_part,weight_mat,p_time)
%instance
%n_square_mat
%P_idx
%job_part
%weight_mat
%p_time

%n_job = instance(1,1);
% instance = instance';
% 
% instance = instance(2:end);

%a = (length(instance)/instance(1,1));

% a_all = {};
% 
% instances  = instance;
% 
% for c = 1:3
%     a1 = instances(1:10); %(1:(length(instance)/instance(1,1))) ; instances(1:(length(instance)/instance(1,1)))
%     instances(1:(length(instance)/instance(1,1))) = [];
%     a_all(end+1) = {a1};
% end

% processing_time = cell2mat(a_all(1,1));
% p_time = processing_time(1,1);
% release_date_all = cell2mat(a_all(1,2))';
% weight_job_all = cell2mat(a_all(1,3))';

%[~,n_job]=size(processing_time);
P_idx = 1:n_job;

%%%%%%p_time = instance(1,1) %Equal Processing Time %%%%%%
n_square_mat = n_job*p_time;
weight_mat = zeros(n_square_mat,n_square_mat);

job_part = repmat(1:p_time, 1, n_job)';

%[~,n_job]=size(processing_time); %size(instance);
%[row,~] = size(processing_time);
row = n_job;
% instance_u_temp = NaN(row,n_job);
% for i = 1:row
%     i;
%     q = instance(i,:);
%     instance_u_temp(i,:) = q(P_idx);
% end

%%%%%%release_date_all = instance(2,:)' %%%%%%

number_job = (1:n_job)';
number_and_release = [number_job release_date_all];

%%dibawah ini sudah bener
%job_part_rd_index = reshape( repmat( [1:n_job], 2,1 ), 1, [] )';
%job_index = reshape( repmat( [1:n_job], 2,1 ), 1, [] )'; %SUDAH BENER

job_part_rd_index = reshape( repmat( P_idx, p_time,1 ), 1, [] )';
job_index = reshape( repmat( P_idx, p_time,1 ), 1, [] )';

for i = 1:n_job
    find(job_part_rd_index == job_index) = number_and_release(job_index,2);
end

release_per_job_part = find';

job_part_release_date = [job_part_rd_index release_per_job_part];

before_order = [job_part_rd_index release_per_job_part job_part];
before_minus = repelem(1,length(before_order))';
before_order_summation = [before_order before_minus];

for i =1:length(before_order_summation)
    x = before_order_summation(i,2) + before_order_summation(i,3) - before_order_summation(i,4);
    before_calculation(i,1) = x;
end

before_summation = [before_order_summation before_calculation];
before_summation = splitvars(table(before_summation));
before_summation.Properties.VariableNames = {'JobIndex' 'JobReleaseTime' 'JobPart' 'Minus1' 'SchedulePlacement'};
fix_before_summation= before_summation;
fix_before_summation.Minus1=[]; %this will remove the "Minus1" column

MatrixIndex = (1:length(before_order))';
%%%%%%weight_job_all = instance(3,:)' %%%%%%

number_and_weight = [number_job weight_job_all];
for i = 1:n_job
    %job_index = number_job(i)
    find(job_part_rd_index == job_index) = number_and_weight(job_index,2);
end

%weight_per_job_part = table(find');
weight_per_job_part = find';
weight_per_job_part = table(weight_per_job_part);
weight_per_job_part = splitvars(table(weight_per_job_part));
weight_per_job_part.Properties.VariableNames = {'JobWeight'};

indexed_fix_before_summation = [table(MatrixIndex) fix_before_summation weight_per_job_part];

for i = 1: length(weight_mat)
    % column 5 = SchedulePlacement
    % column 6 = JobWeight
    % JobWeight multiplied by SchedulePlacement as it was the first time will be placed on time interval
    weight_mat(i,table2array(indexed_fix_before_summation(i,5))) = table2array(indexed_fix_before_summation(i,6)) * table2array(indexed_fix_before_summation(i,5));
end
release_weight_mat = weight_mat;

[~,maxjobparts] = max(indexed_fix_before_summation.JobPart);
[~,minjobparts] = min(indexed_fix_before_summation.JobPart);
%Bottom and upper parts of the jobs
%Loop ke kanan index hasil scheduleplacement + 1
%Loop ke bawah find index for job parts
for i = 1:length(weight_mat)
    indexmatrix = i;
    if table2array(indexed_fix_before_summation(i,4)) == maxjobparts
        rowmatrix = table2array(indexed_fix_before_summation(i,5));
        for w = 1:length(weight_mat)
            w;
            weight_mat(i, rowmatrix + w) = table2array(indexed_fix_before_summation(i,6)) * (rowmatrix + w);
        end
    elseif table2array(indexed_fix_before_summation(i,4)) == minjobparts
        rowmatrix = table2array(indexed_fix_before_summation(i,5));
        for w = 1:length(weight_mat)
            w;
            weight_mat(i, rowmatrix + w) = table2array(indexed_fix_before_summation(i,6)) * (rowmatrix + w);
        end
    end
end
full_weight_mat = weight_mat(1:n_square_mat, 1:n_square_mat);

%Matrix for job assignment at time interval
[rows, ~] = size(weight_mat);
for i = 1:rows
    %Check the column, if full of zero then put "1" on the
    %column, else check next column
    if table2array(indexed_fix_before_summation(i,4)) == maxjobparts
        rowmatrix = table2array(indexed_fix_before_summation(i,5));
        for w = 1:rows
            w;
            assign_mat(i, rowmatrix) = 2;
        end
    elseif table2array(indexed_fix_before_summation(i,4)) == minjobparts
        rowmatrix = table2array(indexed_fix_before_summation(i,5));
        for w = 1:rows
            w;
            assign_mat(i, rowmatrix) = 1;
        end
    end
end
end
