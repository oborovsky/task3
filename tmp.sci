function x = Gauss(A,B)
    C = [A B];
    [n,m] = size(C);
    for i=1:n-1
        Cmax = C(i,i);
        Imax = i;
        for k = i+1:n
            if abs(Cmax) < abs(C(k,i)) then
                Cmax = C(k,i);
                Imax = k;
            end
        end
        if Imax <> i then
            tmp = C(Imax,:);
            C(Imax,:) = C(i,:);
            C(i,:) = tmp;
        end
//        disp(Imax);
//        disp(Cmax);
//        disp(C);
        k = C(i,i);
        if k <> 0 then
            C(i,:) = C(i,:)./k;
        end
        
        a = C(i,:);
        
        for k = i+1:n
            b = C(k,i);
            ak = C(k,:);
            C(k,:) =  ak - a.*b;
        end
        for k = n:-1:2
            b = C(k,k);
            if b <> 0 then
                C(k,:) = C(k,:)./b;
            end
            
            ak = C(k,:);
            for l = k-1:-1:1
                al = C(l,:);
                c = C(l,k);
                C(l,:) = al - ak.*c;
            end
        end
    end
    x = C(:,m)';
endfunction
