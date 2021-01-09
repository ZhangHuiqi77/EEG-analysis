function T = cat_in_holl(n)
    T = zeros(1,n);
    for k = 1:n
        c = unidrnd(3,1);
        while c~= 1
            if c ==2
                T(k) = T(k) + 4;
            else
                T(k) = T(k) + 6;
            end
            c = unidrnd(3,1);
        end
        T(k) = T(k) +2;
    end
end