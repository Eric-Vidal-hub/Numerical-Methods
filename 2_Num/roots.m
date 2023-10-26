clear all;  % Clear all variables
clc;        % Clear the command window

% 1. Create the anonymous function

% An anonymous function is:
% - not stored in a program file (no file name = ”anonymous”);
% - identified by a variable whose data type = ”function handle”;
% - accepting multiple inputs but returning ONE output;
% - allowed to contain only a SINGLE executable statement;
% - inheriting from the current scope any variable
% not in the argument list.
funct = @(x,c) exp(-x) * (x .^ 2+ 2 * x + 2) - c; % funct(x,c) has 2 variables
c=0;
func = @(x) funct(x,c); % func(x) has 1 variable when c is fixed

% 1. Plot the function
for i = 1:21
    x(i) = i - 1;
    y(i) = func(x(i));
end
plot(x,y);      % plot the function

%2. Visualize the function for different c
c = [-1, 0, 1, 1E-6, 0.1, 0.5, 1, 1.99, 2, 2.01, 3]
for k = c
    for i = 1:21
        x(i) = i - 1;
        y(i) = func(x(i));
    end
    plot(x,y);      % plot the function
    hold on;
end
hold off;

% % X. Find the root of the function

% % The function fzero finds a zero of a function.
% % The function fzero requires two inputs:
% % - the function handle of the function to be solved;
% % - an initial guess for the solution.
% % The function fzero returns one output:
% % - the solution of the function.

% x0 = 1; % initial guess
% x = fzero(func,x0); % x is the solution of the function

% % X. Plot the function and the root

% hold on; % hold the current plot
% plot(x,func(x),'ro'); % plot the root
% hold off; % release the current plot

% % X. Find the root of the function for different values of c

% % The function fzero finds a zero of a function.
% % The function fzero requires two inputs:
% % - the function handle of the function to be solved;
% % - an initial guess for the solution.
% % The function fzero returns one output:
% % - the solution of the function.

% c = 0:0.5:5; % c is a vector of 11 points

% for i=1:length(c)
%     func = @(x) funct(x,c(i)); % func(x) has 1 variable when c is fixed
%     x0 = 1; % initial guess
%     x = fzero(func,x0); % x is the solution of the function
%     hold on; % hold the current plot
%     plot(x,func(x),'ro'); % plot the root
%     hold off; % release the current plot
% end