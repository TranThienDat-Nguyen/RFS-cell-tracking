function [w_new,x_new,P_new]= gaus_prune(w,x,P,threshold)

idx= find( w > threshold );
if length(idx) ~= 0
    w_new= w(idx);
    x_new= x(:,idx);
    P_new= P(:,:,idx);
else
    w_new= w;
    x_new= x;
    P_new= P;
end