function [a0, a1, r2] = linreg(x, y)
    % [a0, a1, r2] = linreg(x, y)
    % Found online at: https://www.solvedlib.com/by-using-matlab-linreg-function-file-given-below,186553
    % Rewritten by: Michael Itzkin - September 7, 2022
    %
    % Performs linear regression on the linear x and y dataset
    %
    % Inputs:
    % -x: linear independent data set
    % -y: linear dependent data set
    %
    % Outputs:
    % a0: Constant in y=a1*X + a0
    % a1: Gradient in y=a1*X + a0
    % r2: COefficient of determination

    % Getting the best regression coefficients
    n = length(x);
    sx = sum(x);
    sy = sum(y);
    sx2 = sum(x.^2);
    sxy = sum(x.*y);
    a1 = (n*sxy - sx*sy) / (n*sx2-sx^2);
    a0 = mean(y) - a1*mean(x);

    % Getting R2 value
    st = sum((y - mean(y)).^2);
    sr = sum((y - a0 - a1*x).^2);
    r2 = (st - sr) / st;

end