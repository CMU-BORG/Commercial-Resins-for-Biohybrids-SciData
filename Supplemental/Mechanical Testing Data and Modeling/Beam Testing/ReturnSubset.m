function f_sub = ReturnSubset(F,X,ind)
    f_ = F(X);
    f_sub = f_(ind,:);
end