# spline function from terry :)
function spline!(x,y,n,yp1,ypn,y2)
    
    u = zeros(n)
    if (yp1 > 0.99e30)
        y2[1] = 0.0
        u[1] = 0.0
    else
        y2[1] = -0.5
        u[1] = (3.0/(x[2] - x[1])) * ( (y[2] - y[1]) / (x[2] - x[1]) - yp1)
    end
    
    for i=2:n-1
        sig = (x[i] - x[i-1])/(x[i+1]-x[i-1])
        p = sig * y2[i-1] + 2.0
        y2[i]=(sig-1.0)/p
        u[i] = (6.0 * ((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1]) /
            (x[i]-x[i-1]))/(x[i+1]-x[i-1]) - sig * u[i-1])/p
    end
    
    if (ypn > 0.99e30)
        qn=0.0
        un=0.0
    else
        qn = 0.5
        un = (3.0/(x[n]-x[n-1])) * (ypn - (y[n]-y[n-1])/(x[n]-x[n-1]))
    end
    
    y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0)
    for k=n-1:-1:1
        y2[k]=y2[k]*y2[k+1]+u[k]
    end

end

function splint(xa,ya,y2a,n,x)
    
    klo = 1
    khi = n
    while (khi-klo > 1)
        k = Int(round((khi + klo)/2))
        if (xa[k] > x) 
            khi=k
        else
            klo=k
        end
    end
    h = xa[khi]-xa[klo]
    if (h == 0.0) 
        println("bad xa input in splint")
    end
    a = (xa[khi]-x)/h
    b = (x - xa[klo])/h
    y = a * ya[klo] + b*ya[khi] + ((a^3.0 - a)*y2a[klo]+(b^3.0-b)*y2a[khi])*(h^2.0)/6.0

    return y

end 
