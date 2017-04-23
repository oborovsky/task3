deff('y=u(x)','y=(1-exp(100*x))/(1-exp(100))');

function x=shufle(m, d)
    a = m(1,:);
    b = m(2,:);
    c = m(3,:);
    s = length(b);
    et(1) = 0;
    ks(1) = 0;
    
    for i=1:s
        ks(i+1) = (-c(i)) / (a(i)*ks(i) + b(i));
        et(i+1) = (d(i) - a(i) * et(i)) / (a(i) * ks(i) + b(i));
    end
    
    x(s+1) = 0;
    
    for j=0:s-1
        x(s-j) = ks(s-j+1)*x(s-j+1) + et(s-j+1);
    end
    x(s+1) = [];
endfunction
// аппроксимация для (5)
function [m, d] = makeLinearSystem1 (s,h)
    a = [0];
    b = [1];
    c = [0];
    dd = [0];
   
    k1 = -(1/h^2 + 100/(2*h));
    k2 = 2/h^2;
    k3 = (100/(2*h) - 1/h^2);
    
    for i=2:s-1
        a(i) = k1;
        b(i) = k2;
        c(i) = k3;    
        dd(i) = 0;
    end

    a(s) = 0;
    b(s) = 1;
    c(s) = 0;
    dd(s) = 1;
    
    m = [a';b';c'];
    d = dd';
endfunction
// аппроксимация для (6)
function [m, d] = makeLinearSystem2 (s,h)
    a = [0];
    b = [1];
    c = [0];
    dd = [0];
  
    k1 = -(1/h^2 + 100/(h));
    k2 = 2/h^2 + 100/h;
    k3 = -1/h^2;
    
    for i=2:s-1
        a(i) = k1;
        b(i) = k2;
        c(i) = k3;    
        dd(i) = 0;
    end

    a(s) = 0;
    b(s) = 1;
    c(s) = 0;
    dd(s) = 1;
    
    m = [a';b';c'];
    d = dd';
endfunction

function [e1,e2,Uex,y1,y2] = makeApp (X)
    h = 1/X;
    x = 0:h:1;
    Uex = u(x);
    s = length(x);
    
    [m1,d1] = makeLinearSystem1(s, h);
    [m2,d2] = makeLinearSystem2(s, h);
    
    tmp = shufle(m1, d1);
    y1 =  tmp';
    tmp = shufle(m2, d2);
    y2 = tmp';
   
    e1 = y1 - Uex;
    e2 = y2 - Uex;
endfunction

N = 1000;
[e1,e2,Uex,y1,y2] = makeApp(N);
ee1(1) = max(abs(e1));
ee2(1) = max(abs(e2));
for i = 2:2
    N = 2*N;
    X = N;
    [e1,e2,Uex,y1,y2] = makeApp(X);
    ee1(i) = max(abs(e1));
    ee2(i) = max(abs(e2));
    x = 0:1/X:1;
    plot(x,Uex,'-b');
    plot(x, y1,'-g');
    plot(x,y2,'-r');
    printf("ee1(%d)=%f,ee1(%d)=%f, p = %f\n",i-1,ee1(i-1),i,ee1(i), abs(log2(ee1(i-1)/ee1(i))));
    printf("ee2(%d)=%f,ee2(%d)=%f, p = %f\n",i-1,ee2(i-1),i,ee2(i), abs(log2(ee2(i-1)/ee2(i))));
end
