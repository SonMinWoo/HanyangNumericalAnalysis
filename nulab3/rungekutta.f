implicit real*8(a-h, o-z)
real*8 k1, k2, k3, k4

exactly(x) = x/(1+x**2)
f(x,y) = 1/(1+x**2)-2*y**2

t0 = 0; b = 1; h = 0.2; t = t0
y0 = 0.0; w0 = y0

write(*,*)
write(*,*) 'When h = 0.2 ------------------------'

do while (t .le. b)
    t = t0 + h
    k1 = h*f(t0, w0)
    k2 = h*f(t0 + h/2, w0 + k1/2)
    k3 = h*f(t0 + h/2, w0 + k2/2)
    k4 = h*f(t, w0 + k3)
    wnew = w0 + (k1 + 2*k2 + 2*k3 + k4)/6 ! RK4 의 근사항
    fval = exactly(t)

    write(*,*) 'Runge-Kutta value and Exact value (x, induced value, exact value) : ', t, wnew, fval

    d = ABS(wnew - fval)

    t0 = t
    w0 = wnew
end do

write(*,*)
write(*,*) '(x, induced value of Runge-Kutta, exact value, difference) : ', t, wnew, fval, d

read(*,*)
stop
