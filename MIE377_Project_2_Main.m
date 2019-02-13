%% MIE377 (Winter 2018) - Project 2
% The purpose of this program is to test the out-of-sample performance of 
% 5 different portfolio optimization models. We will test the following
% models:
% 1. Nominal MVO
% 2. Robust MVO
% 3. Resampling MVO
% 4. Most Diverse MVO
% 5. Conditional Value at Risk
%
%
% Student Name: Shen Lin
% Student ID: 1001367193
% Student Name: Zhuoran Lu
% Student ID: 1002129127
% Student Name: Cong Shi
% Student ID: 1002330587

clc
clear all
format long
% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that Matlab 2017b version doesn't support datetime format. The date format of excel
% has been changed to "long-date". E.g 1/6/2012 has been changed to 6-Jan-12


% Load the stock weekly prices and factors weekly returns
adjClose = readtable('Project1_Data_adjClose.csv', 'ReadRowNames', true);
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Properties.RowNames));

factorRet = readtable('Project1_Data_FF_factors.csv', 'ReadRowNames', true);
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of shares outstanding per company
mktNoShares = [36.40 9.86 15.05 10.91 45.44 15.43 43.45 33.26 51.57 45.25 ...
                69.38 8.58 64.56 30.80 36.92 7.70 3.32 54.68 38.99 5.78]' * 1e8;

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Investment strategies
% Note: You must populate the functios MVO.m, MVO_card.m and BL.m with your
% own code to construct the optimal portfolios. 
funNames  = {'MVO' 'Robust MVO' 'Resampling MVO' 'Most Diverse MVO' 'CVaR'};
NoMethods = length(funNames);

funList = {'MVO' 'Robust_MVO' 'Resampling_MVO' 'Diverse_MVO','cvar'};
funList = cellfun(@str2func, funList, 'UniformOutput', false);

% Define parameters for robust MVO
beta_robust = 0.9;

% Define parameters for resampling MVO
T = 100;
N = 60;

% Define parameters for most diverse MVO
n_assets = 12;

% Define parameters for CVaR
no_scenarios = 2000;
beta_cvar = 0.95;
dt = 26;

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns, cov. matrix,
% etc) from the Fama-French factor models. You will have to re-estimate 
% your parameters at the start of each rebalance period, and then 
% re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;
currentVal = NaN(NoPeriods,NoMethods);
currentVal1 = NaN(NoPeriods,5);
currentVal2 = NaN(NoPeriods,5);
currentVal3 = NaN(NoPeriods,5);
currentVal4 = NaN(NoPeriods,5);
currentVal5 = NaN(NoPeriods,5);

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        
        currentVal(t,:) = initialVal;
        
    else
        for i = 1 : NoMethods
            
            currentVal(t,i) = currentPrices' * NoShares{i};
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % Calculate your initial exp. returns and covariance matrix using the
    % Fama-French factor model. You may find these values as you prefer,
    % the code given below is just an example. 
    n=size(tickers, 1);
    Beta_co=zeros(4, n);%four coefficients: a, beta_im, beta_is, beta_iv
    num_weeks = size(periodReturns, 1);
    for i = 1:n
        Beta_co(:, i) = regress(periodReturns(:,i),[ones(num_weeks,1),periodFactRet]); 
    end
    Alpha = Beta_co(1,:)';          % n x 1 vector of alphas
    V = Beta_co(2:end, :);          % m x n matrix of betas
    f_bar = (geomean(1+ periodFactRet)-1)';      % m x 1 vector of factor expected returns
    F = cov(periodFactRet);          % m x m factor covariance matrix
    epsilon = periodReturns - [ones(num_weeks, 1) periodFactRet] * Beta_co;    % Regression residuals
    mu = (geomean(1 + periodReturns, 1) - 1)';    % n x 1 vector of asset exp. returns
    Q  = cov(periodReturns);         % n x n asset covariance matrix
    D = cov(epsilon);          % Diagonal n x n matrix of residual variance
    %----------------------------------------------------------------------
    
    % Define the target return for the 2 MVO portfolios
    targetRet = mean(mu);
    
    % Optimize your portfolios to get the weights 'x'
    x{1}(:,t) = funList{1}(mu, Q, targetRet); 
    x{2}(:,t) = funList{2}(periodReturns, Q, beta_robust);  
    x{3}(:,t) = funList{3}(mu, Q, targetRet, T, N); 
    x{4}(:, t) = funList{4}(mu, Q, targetRet, n_assets); 
    [x{5}(:, t), Ret] = funList{5}(mu, Q, targetRet, no_scenarios, currentPrices, beta_cvar, dt); 
    
    MVO_cvar(t) = model_cvar(x{1}(:,t), Ret, 0.95);
    Robust_cvar(t) = model_cvar(x{2}(:,t), Ret, 0.95);
    Resampling_cvar(t) = model_cvar(x{3}(:,t), Ret, 0.95);
    Diverse_cvar(t) = model_cvar(x{4}(:,t), Ret, 0.95);
    cvar_cvar(t) = model_cvar(x{5}(:,t), Ret, 0.95);
    
    
    % Calculate the optimal number of shares of each stock you should hold
    
    for i = 1 : NoMethods
        
        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t,i) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,i) = periodPrices * NoShares{i};
        
        % *************** WRITE YOUR CODE HERE ***************
        %------------------------------------------------------------------
        
        % Calculate your transaction costs for the current rebalance
        % period. The first period does not have any cost since you are
        % constructing the portfolios for the 1st time. 
        
        if t ~= 1
          
           tCost(t-1, i) = 0.005 * abs(NoSharesOld{i} - NoShares{i})' * currentPrices;
        end
        
        NoSharesOld{i} = NoShares{i};
        %------------------------------------------------------------------
        
    end
    
    % Calculate the Ex Ante Sharp Ratio    
    for i = 1:size(funList,2)

    SharpeRatio_ExAnte(i,t) = mu'*x{i}(:,t) / sqrt(x{i}(:,t)' * Q * x{i}(:,t));
    
    end
    
    % Calculate the Ex Post Sharp Ratio   
    ret_p = ( portfValue(fromDay+1:toDay,:) - portfValue(fromDay:toDay-1,:) )...
    ./ portfValue(fromDay:toDay-1,:);
    
    risk_free_p = riskFree((testStart <= dates & dates <= testEnd), :);
    risk_free_p = risk_free_p(2:end,:);
    ret_p_excess = ret_p - (diag( table2array(risk_free_p)) * ones( size(risk_free_p)));

    std_p_excess = sqrt(var(ret_p_excess));
    mu_p_excess = geomean(1+ret_p_excess) - 1;
    SharpeRatio_ExPost(:,t) = mu_p_excess ./ std_p_excess;
    

    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);
   
    
    % study the effect of having different confidence levels for robust MVO
    beta = [0.75, 0.8, 0.85, 0.9, 0.95];
    if t == 1  
        currentVal1(t,:) = initialVal;
    else
        for i = 1 : length(beta)    
            currentVal1(t,i) = currentPrices' * NoShares_robust{i};
        end
    end
        
    for i=1:length(beta)        
        x_robust{i}(:,t) = funList{2}(periodReturns, Q, beta(i));
        NoShares_robust{i} = x_robust{i}(:,t) .* currentVal1(t,i) ./ currentPrices;
        portfValue_robust(fromDay:toDay,i) = periodPrices * NoShares_robust{i};
    end
    
    % study the effect of having different confidence levels for CVaR
    if t == 1  
        currentVal2(t,:) = initialVal;
    else
        for i = 1 : length(beta)    
            currentVal2(t,i) = currentPrices' * NoShares_cvar{i};
        end
    end
    
    for i=1:length(beta)
        x_cvar{i}(:,t) = funList{5}(mu, Q, targetRet, 2000, currentPrices,beta(i), dt);
 %        x_cvar{i}(:,t) =  funList{1}(mu, Q, targetRet); 
        NoShares_cvar{i} = x_cvar{i}(:,t) .* currentVal2(t,i) ./ currentPrices;
        portfValue_cvar(fromDay:toDay,i) = periodPrices * NoShares_cvar{i};
    end
    
    % study the effect of having different number of simulations for
    % resampling
    no_simulations_resamp = [20, 60, 100, 500, 1000];
    if t == 1  
        currentVal3(t,:) = initialVal;
    else
        for i = 1 : length(no_simulations_resamp)    
            currentVal3(t,i) = currentPrices' * NoShares_resamp{i};
        end
    end
    
    for i=1:length(no_simulations_resamp)
        x_resamp{i}(:,t) = funList{3}(mu, Q, targetRet, no_simulations_resamp(i), 60);
        NoShares_resamp{i} = x_resamp{i}(:,t) .* currentVal3(t,i) ./ currentPrices;
        portfValue_resamp(fromDay:toDay,i) = periodPrices * NoShares_resamp{i};
    end
    
    % study the effect of having different number of simulations for
    % CVaR
    no_simulations_cvar = [500, 1000, 1500, 2000, 2500];
    if t == 1  
        currentVal4(t,:) = initialVal;
    else
        for i = 1 : length(no_simulations_resamp)    
            currentVal4(t,i) = currentPrices' * NoShares_cvar_no{i};
        end
    end
    
    for i=1:length(no_simulations_cvar)
        x_cvar_no{i}(:,t) = funList{5}(mu, Q, targetRet, no_simulations_resamp(i), currentPrices, 0.95, 26);
        NoShares_cvar_no{i} = x_cvar_no{i}(:,t) .* currentVal4(t,i) ./ currentPrices;
        portfValue_cvar_no(fromDay:toDay,i) = periodPrices * NoShares_cvar_no{i};
    end
    
    % study the effect of having different step size (dt) for CVaR
    dt1 = [1, 6.5, 13, 26];
    if t == 1  
        currentVal5(t,:) = initialVal;
    else
        for i = 1 : length(dt1)    
            currentVal5(t,i) = currentPrices' * NoShares_cvar_dt{i};
        end
    end
    
    for i=1:length(dt1)
        x_cvar_dt{i}(:,t) = funList{5}(mu, Q, targetRet, no_scenarios, currentPrices, beta_cvar, dt1(i));
        NoShares_cvar_dt{i} = x_cvar_dt{i}(:,t) .* currentVal5(t,i) ./ currentPrices;
        portfValue_cvar_dt(fromDay:toDay,i) = periodPrices * NoShares_cvar_dt{i};
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************** WRITE YOUR CODE HERE ***************
%--------------------------------------------------------------------------

% Calculate the portfolio average return, variance (or standard deviation),
% or any other performance and/or risk metric you wish to include in your
% report.

plotDates = dates(datetime('2013-01-01') <= dates);
% portfolio return
portfRet = ( portfValue(2:end,:) - portfValue(1:end-1,:) ) ./ portfValue(1:end-1,:);

% portfolio expected return
mu_p = geomean(portfRet + 1,1) - 1;
% portfolio variances
var_p = cov(portfRet);
% portfolio standard deviation 
sigma_p = std(portfRet);


%% 

%--------------------------------------------------------------------------
% 4.1 Plot the portfolio values 
%--------------------------------------------------------------------------

fig1 = figure(1);
plot(plotDates, portfValue(:,1))
hold on
plot(plotDates, portfValue(:,2))
hold on
plot(plotDates, portfValue(:,3))
hold on
plot(plotDates, portfValue(:,4))
hold on
plot(plotDates, portfValue(:,5))
legend(funNames, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio Value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'Portfolio_Value','-dpng','-r0');




%--------------------------------------------------------------------------
% 4.2 Plot the portfolio weights 
%--------------------------------------------------------------------------

% MVO Plot
fig2 = figure(2);
area(x{1}', 'FaceColor','flat');
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .png for use in MS Word
print(fig2,'CVAR for models','-dpng','-r0');


% Robust MVO Plot
fig3 = figure(3);
area(x{2}', 'FaceColor','flat');
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Robust MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig3,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .png for use in MS Word
print(fig3,'MVO_Robust Weights','-dpng','-r0');

% Resampling MVO Plot
fig4 = figure(4);
area(x{3}', 'FaceColor','flat');
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Resampling MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig4,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig4,'Position');
set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .png for use in MS Word
print(fig4,'MVO_Resample Weights','-dpng','-r0');


% Most_Diverse MVO Plot
fig5 = figure(5);
area(x{4}', 'FaceColor','flat');
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Most Diverse MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig5,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig4,'Position');
set(fig5,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .png for use in MS Word
print(fig5,'MVO_MostDiverse Weights','-dpng','-r0');

% Cvar Plot
fig6 = figure(6);
area(x{5}', 'FaceColor','flat');
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('CVaR Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig6,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig4,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .png for use in MS Word
print(fig6,'Cvar_Weights','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.3 Plot the transaction cost values 
%--------------------------------------------------------------------------

% Plot the transaction costs
fig7 = figure(7);
plot([1:5], tCost(:,1), '-o');
hold on
plot([1:5], tCost(:,2),'-o');
hold on
plot([1:5], tCost(:,3), '-o');
hold on
plot([1:5], tCost(:,4), '-o');
hold on
plot([1:5], tCost(:,5), '-o');
legend(funNames, 'Location', 'eastoutside','FontSize',12);
title('Transaction Cost', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);
xlabel('n-th rebalancing portfolio','interpreter','latex','FontSize',12)

% If you want to save the figure as .png for use in MS Word
print(fig7,'Transaction Cost','-dpng','-r0');
%% 

%--------------------------------------------------------------------------
% 4.4 Plot the effect of having different parameters 
%--------------------------------------------------------------------------

% Plot effect of having different confidence levels for robust MVO
fig8 = figure(8);
for i = 1:size(portfValue_robust, 2)
    plot(plotDates, portfValue_robust(:,i))
    hold on
end
leg = {'\beta = 0.75', '\beta = 0.8', '\beta = 0.85', '\beta = 0.9', '\beta = 0.95'};
legend(leg, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Robust MVO Port. Val. with different beta', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

print(fig8,'robust MVO confidence','-dpng','-r0');

% Plot effect of having different confidence levels for CVaR 
fig9 = figure(9);
for i = 1:size(portfValue_cvar, 2)
    plot(plotDates, portfValue_cvar(:,i))
    hold on
end
legend(leg, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('CVaR Port. Val. with different beta', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

print(fig9,'cvar confidence','-dpng','-r0');
%% 

% Plot effect of having different number of simulations for Resampling
fig10 = figure(10);
for i = 1:size(portfValue_resamp, 2)
    plot(plotDates, portfValue_resamp(:,i))
    hold on
end
leg2 = {'# Scenarios = 20', '# Scenarios = 60', '# Scenarios = 100', '# Scenarios = 500', '# Scenarios = 1000'};
legend(leg2, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Resampling Port, Val. with different # of simulations', 'FontSize', 10)
ylabel('Value','interpreter','latex','FontSize',12);

fig10.PaperUnits = 'inches';
fig10.PaperPosition = [0 0 8 6];
print(fig10,'resamp no sim','-dpng','-r0');

% Plot effect of having different number of simulations for CVaR
fig11 = figure(11);
for i = 1:size(portfValue_cvar_no, 2)
    plot(plotDates, portfValue_cvar_no(:,i))
    hold on
end
leg3 = {'# Scenarios = 500', '# Scenarios = 1000', '# Scenarios = 1500', '# Scenarios = 2000', '# Scenarios = 2500'};
legend(leg3, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('CVaR Port, Val. with different # of simulations', 'FontSize', 10)
ylabel('Value','interpreter','latex','FontSize',12);

print(fig11,'cvar no sim','-dpng','-r0');

% Plot effect of having different number of simulations for CVaR
fig19 = figure(19);
for i = 1:size(portfValue_cvar_dt, 2)
    plot(plotDates, portfValue_cvar_dt(:,i))
    hold on
end
leg3 = {'dt = 1', 'dt = 6.5', 'dt = 13', 'dt = 26'};
legend(leg3, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('CVaR Port, Val. with different # of simulations', 'FontSize', 10)
ylabel('Value','interpreter','latex','FontSize',12);

print(fig19,'cvar dt','-dpng','-r0');

%% 

%--------------------------------------------------------------------------
% 4.5 Plot Sharpe Ratios
%--------------------------------------------------------------------------
% Comparison of Ex Ante and Ex Post for each portfolio
for i = 1:NoMethods
    fig12{i} = figure(i+11);
    plot([1:6], SharpeRatio_ExAnte(i, :))
    hold on
    plot([1:6], SharpeRatio_ExPost(i, :))
    
    legend({'Ex Ante', 'Ex Post'}, 'Location', 'eastoutside','FontSize',12);
%     datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
    set(gca,'XTickLabelRotation',30);
    tit = strcat(funNames{i}, ' Sharpe Ratio');
    title(tit, 'FontSize', 10)
    ylabel('SR','interpreter','latex','FontSize',12);

%     Define the plot size in inches
    set(fig12{i},'Units','Inches', 'Position', [0 0 8, 5]);
    pos1 = get(fig12{i},'Position');
    set(fig12{i},'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos1(3), pos1(4)]);

    print(fig12{i}, tit,'-dpng','-r0');
end
%% 

% Comparison of Ex Ante
for i = 1:NoMethods
    fig17 = figure(17);
    plot([1:6], SharpeRatio_ExAnte(i, :))
    hold on
end

legend(funNames, 'Location', 'eastoutside','FontSize',12);
%     datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Sharpe Ratio Ex Ante', 'FontSize', 10)
ylabel('SR','interpreter','latex','FontSize',12);

%     Define the plot size in inches
set(fig17,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig17,'Position');
set(fig17,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos1(3), pos1(4)]);
print(fig17, 'Sharpe Ratio Ex Ante', '-dpng','-r0');

% Comparison of Post Ante
for i = 1:NoMethods
    fig18 = figure(18);
    plot([1:6], SharpeRatio_ExPost(i, :))
    hold on
end

legend(funNames, 'Location', 'eastoutside','FontSize',12);
%     datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Sharpe Ratio Ex Post', 'FontSize', 10)
ylabel('SR','interpreter','latex','FontSize',12);

%     Define the plot size in inches
set(fig18,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig18,'Position');
set(fig18,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos1(3), pos1(4)]);
print(fig18, 'Sharpe Ratio Ex Post', '-dpng','-r0');

%--------------------------------------------------------------------------
% 4.6 Plot the CVaR values for each model
%--------------------------------------------------------------------------
fig20 = figure(20);
plot([1:6], MVO_cvar, '--');
hold on
plot([1:6], Robust_cvar,'--');
hold on
plot([1:6], Resampling_cvar, '--');
hold on
plot([1:6], Diverse_cvar, '--');
hold on
plot([1:6], cvar_cvar, '--');
legend(funNames, 'Location', 'eastoutside','FontSize',12);
title('CVaR Value for each model', 'FontSize', 14)
ylabel('CVaR value','interpreter','latex','FontSize',12);
xlabel('Rebalance Periods','interpreter','latex','FontSize',12)

% If you want to save the figure as .png for use in MS Word
print(fig20,'CVaR_value','-dpng','-r0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End