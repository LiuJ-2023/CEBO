function f_best = CEBO_CEC2010(problem_num)
if problem_num == 1
    up = [10 10 10 10 10 10 10 10 10 10];
    dn = [0 0 0 0 0 0 0 0 0 0];
    d = 10;
    df = 3;
elseif problem_num == 7
    up = 140*ones([1,10]);
    dn = -140*ones([1,10]);
    d = 10;
    df = 2;
elseif problem_num == 8
    up = 140*ones([1,10]);
    dn = -140*ones([1,10]);
    d = 10;
    df = 2;
elseif problem_num == 13
    up = 500*ones([1,10]);
    dn = -500*ones([1,10]);
    d = 10;
    df = 4;
elseif problem_num == 14
    up = 1000*ones([1,10]);
    dn = -1000*ones([1,10]);
    d = 10;
    df = 4;
elseif problem_num == 15
    up = 1000*ones([1,10]);
    dn = -1000*ones([1,10]);
    d = 10;
    df = 4;
end
% Parameter Setting
pop_size = 30;
alpha = 0.05;
pop = (up-dn).*rand([pop_size,d]) + dn;

% Initialize Population and Database
ND = 50;
DB_x = (up-dn).*rand([ND,d]) + dn;
DB_y = Fitness_2010(DB_x,problem_num);

% Recored the best solution
[x_best,f_best] = Selection_FR(DB_x,DB_y);

% Main loop
for gen = 1:25  
    % Build GP Surrogates for the objective functions and constraints
    for j = 1:df
        dmodel{j} = fitrgp(DB_x,DB_y(:,j),'KernelFunction','squaredexponential','Sigma',0.01);
    end
    
    %% First Strategy
    % Search for the query solution
    x_q = DE_FR_opt(dmodel,up,dn,alpha,pop);
    f_q = Fitness_2010(x_q,problem_num);
    % Update the Database
    DB_x = [DB_x;x_q];
    DB_y = [DB_y;f_q];
    % Update alpha
    CV_q = sum(max(f_q(2:end),0));
    if CV_q <= 0
        alpha = 0.9*alpha;
    end
    % Update the current best solution
    [x_best,f_best] = Selection_FR(DB_x,DB_y);
    % Print
    clc
    fprintf(['当前最优解目标函数值为',num2str(f_best(1)),'\n'])
    fprintf(['当前最优解约束违反量为',num2str(sum(max(f_best(2:end),0))),'\n'])
    fprintf(['当前插入解目标函数值为',num2str(f_q(1)),'\n'])
    fprintf(['当前插入解约束违反量为',num2str(sum(max(f_q(2:end),0))),'\n'])
    fprintf(['内点参数',num2str(alpha),'\n'])
    
    %% Second Strategy
    % Evolution of the population
    for i = 1:pop_size
        % Generate offspring
        u = crossover(pop,up,dn,10,i);
        x = [pop(i,:);u];
        for j = 1:df
            if j == 1
                [f_hat_x(:,j),s_hat_x(:,j)] = predict(dmodel{j}, x); 
                f_cb(:,j) = f_hat_x(:,j) - 2*s_hat_x(:,j);
            else
                [f_hat_x(:,j),s_hat_x(:,j)] = predict(dmodel{j}, x); 
                f_cb(:,j) = f_hat_x(:,j) - 2*s_hat_x(:,j);
            end
        end
        % Update population by using feasibility rule
        [off(i,:),f_off(i,:)] = Selection_FR(x,f_cb);
    end
    pop = off;
    f = f_off;
    % Randomly select a solution from the population  
    [len_popest,~] = size(pop);
    idx_sel = ceil(len_popest*rand);
    x_q = pop(idx_sel,:);
    f_q = Fitness_2010(x_q,problem_num);
    % Update the database
    DB_x = [DB_x;x_q];
    DB_y = [DB_y;f_q];
    % Update the current best solution
    [x_best,f_best] = Selection_FR(DB_x,DB_y);   
    
    % Print
    fprintf(['当前迭代数',num2str(gen),'\n'])
    fprintf(['当前最优解目标函数值为',num2str(f_best(1)),'\n'])
    fprintf(['当前最优解约束违反量为',num2str(sum(max(f_best(2:end),0))),'\n'])
    fprintf(['当前插入解目标函数值为',num2str(f_q(1)),'\n'])
    fprintf(['当前插入解约束违反量为',num2str(sum(max(f_q(2:end),0))),'\n'])
    fprintf(['内点参数',num2str(alpha),'\n'])
    
end
 
function offspring = crossover(x,up,dn,total_num,i)
[m,d] = size(x);
for num = 1:total_num
    % Uniform crossover
    mask = (rand(1,d) <= 0.5);
    rs = randperm(m,1);
    offspring(num,:) = x(i,:).*mask + x(rs,:).*(1-mask);
    % Gaussian mutation
    offspring(num,:) = offspring(num,:) + 0.01*(up - dn).*randn(1,d);
    % Repair
    offspring(num,:) = (offspring(num,:) >= up).*max(dn,2*up-offspring(num,:)) + (offspring(num,:) < up).*offspring(num,:);
    offspring(num,:) = (offspring(num,:) <= dn).*min(up,2*dn-offspring(num,:)) + (offspring(num,:) > dn).*offspring(num,:);
end

function u = DE(x,up,dn,total_num,i)
[m,d] = size(x);
for num = 1:total_num
    rs = randperm(m,3);
    rj = rand(1,d);
    % DE/rand/1
    v(num,:) = x(rs(1),:) + 0.5*(x(rs(2),:) - x(rs(3),:));
    % bin交叉
    u(num,:) = v(num,:).*(rj<0.9) + x(i,:).*(rj>=0.9);
    % 修补
    u(num,:) = (u(num,:) >= up).*max(dn,2*up-u(num,:)) + (u(num,:) < up).*u(num,:);
    u(num,:) = (u(num,:) <= dn).*min(up,2*dn-u(num,:)) + (u(num,:) > dn).*u(num,:);
end

function [x_best,f_best] = Selection_FR(x,f_x)
F = f_x(:,1);
CV = sum(max(f_x(:,2:end),0),2);

idx_fea = (CV <= 0);
if sum(idx_fea) == 0
    [~,idx_best] = min(CV);
    x_best = x(idx_best,:);
    f_best = f_x(idx_best,:);
else
    F_FR = idx_fea.*F + ~idx_fea.*( max(F(idx_fea)) + CV );
    [~,idx_best] = min(F_FR);
    x_best = x(idx_best,:);
    f_best = f_x(idx_best,:);
end

function [x_best] = DE_FR_opt(dmodel,up,dn,g_rep,initpop)
d = length(up);
pop1 = initpop;
pop2 = (up - dn).*rand([25,d]) + dn;
[sizee,~] = size(pop1);
idx_selperm = randperm(sizee,sizee/2);
pop = pop1(idx_selperm,:);
pop = [pop1;pop2];
[pop_size,~] = size(pop);

for i = 1:pop_size
    for j = 1:length(dmodel)
        if j == 1
            [v(i,j),s(i,j)] = predict(dmodel{j}, pop(i,:));
            f(i,j) = v(i,j) - 2*s(i,j);
        else
            [v(i,j),s(i,j)] = predict(dmodel{j}, pop(i,:));
            f(i,j) = v(i,j) + g_rep*s(i,j);
        end
    end
end

for gen = 1:100
    for i = 1:pop_size
        [x_best,f_best] = Selection_FR(pop,f);
        u(i,:) = DE(pop,up,dn,1,i);
        for j = 1:length(dmodel)
            if j == 1
                [v_u(i,j),s_u(i,j)] = predict(dmodel{j}, u(i,:));
                f_u(i,j) = v_u(i,j) - 2*s_u(i,j);
            else
                [v_u(i,j),s_u(i,j)] = predict(dmodel{j}, u(i,:));
                f_u(i,j) = v_u(i,j) + g_rep*s_u(i,j);
            end 
        end
        x = [pop(i,:);u(i,:)];
        f_x = [f(i,:);f_u(i,:)];
        [pop(i,:),f(i,:)] = Selection_FR(x,f_x);
    end
end
[x_best,f_best] = Selection_FR(pop,f);
