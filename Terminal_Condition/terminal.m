function [time,portion] = terminal(x,tol)
if x(end,2)>tol 
    time=x(end,1);
else
    for i=size(x,1):-1:2
        if x(i,2)<=tol && x(i-1,2)>=tol;
            time=x(i,1);
            break;
        end
    end
    if i==2 
        time=0;
    end
end
portion=0;
for i=1:size(x,1)
    if x(i,2)>=tol
        portion=portion+0.01;
    end
end
    