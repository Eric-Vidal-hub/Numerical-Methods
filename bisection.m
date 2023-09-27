function [approx,iter] = bisection(a,b,f,accuracy)
%% BISECTION ALGORITHM 
%  for searching a simple root of f(x) in interval [a,b]
%  such that sign(f(a)) = -sign(f(b)) 
%  while minimizing the number of function evaluations.
%  INPUT:
%  a = left  boundary of search interval
%  b = right boundary of search interval
%  f = function featuring a simple root in [a,b]
%  accuracy = absolute error on root position 
%             (algorithm changes sign if entered negative)
%  OUTPUT:
%  approx = approximation of simple root
%         = NaN if sign(f(a)) = sign(f(b))
%  iter   = number of bisections; 
%         = 0 if sign(f(a)) = sign(f(b))
%  iter+2 = number of function evaluations if sign(f(a))=-sign(f(b)) 
iter=0;  % Iteration counter
fa=f(a); % Avoid recomputing f(x) in two next tests
if( sign(fa)*sign(f(b)) >= 0 ) 
  fprintf('Root is not surely bracketed by input [a,b]\n');
  approx=0/0; % Return a meaningless output
  return
end	  
if(fa < 0)    % Orienting the search so that f(x)>0 lies at x+delta 
  approx= a;  % in while loop below
  delta = b-a;
else
  approx= b;
  delta = a-b;
end
while (abs(delta) > abs(accuracy) ) % Bisection loop
  delta = delta/2; 
  x = approx+delta;
  iter  = iter+1; 
  fbisec = f(x);   % Avoid recomputing f(x) in two next tests
  if(fbisec <= 0)
    approx = x;
  end
  if(fbisec == 0) 
      fprintf('Exact root within machine precision!\n');
    break  
  end
end
end % End of function bisection