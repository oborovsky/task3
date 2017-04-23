deff('y=u(x)','y=(1-exp(100*x))/(1-exp(100))');
deff('y=f0(x)','y=x');
//deff('y=fk(x,k)','y=(exp(k*x)-exp(k))/(1-exp(k)) - (1-x)');
//deff('y=dfk(x,k)','y=k*exp(k*x)/(1-exp(k)) + 1');
//deff('y=d2fk(x,k)','y=k^2*exp(k*x)/(1-exp(k))');
deff('y=Lfk(x,k)','y=-d2fk(x,k) + 100*dfk(x,k)');

deff('y=fk(x,k)','y=1 - cos(2*%pi*k*x)');
deff('y=dfk(x,k)','y=2*%pi*k*sin(2*%pi*k*x)');
deff('y=d2fk(x,k)','y=(2*%pi*k)^2*cos(2*%pi*k*x)');

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

function y = yk(x,c)
    y = f0(x);
    for k=1:length(c)
        y = y + c(k)*fk(x,k)
    end
endfunction

N = 2;
ll = 5;
for l=1:ll
    N = 2*N;
    h = 1/N;
    xn = h:h:1-h;
    A = zeros(N-1,N-1);
    
    for i = 1:N-1
        for j = 1:N-1
            A(i,j) = Lfk(xn(i),j);
        end
    end
    for i = 1:N-1
        F(i) = -100;
    end

    c = Gauss(A,F);
    
    x = 0:0.01:1;
    Uex = u(x);
    y1(l,:) = yk(x,c);

    ee(l) = max(abs(Uex-y1(l,:)));
    NN(l) = N; 
//    printf("e = %f \n",max(abs(Uex-y1)));   
    
end

for l =1:ll
    plot(x,y1(l,:),'-r');
end
plot(x,Uex,'-b');
