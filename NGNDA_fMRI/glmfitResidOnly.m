function resid = glmfitResidOnly(X,y)
% Purpose: Run general linear model only providing residuals. Makes it possible to run par
%          for loop and speed up data processing
% Inputs: X = design matrix of motion regressors
%         y = response variable
% Outputs: resid = residual of data after regressing out motion or whatever
%                  parameters. Desired signal
[~,~,stats] = glmfit(X,y);
resid = stats.resid;
end