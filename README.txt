This is an offline, desktop PyQt version of the online expression evaluator at https://garland.pythonanywhere.com/.  
The offline version offers all the features of the online version as well as animated plots and keyboard shortcuts. It can evaluate and simplify various numerical and algebraic expressions, including substitutions, 
solve polynomial equations and inequalities,  solve systems of linear equations, evaluate algebraic differentiation and numerical integration, handle vector calculus, discrete Fourier transform, numerically solve 
ordinary differential equations, including boundary value problems
If you have Python3 installed, you can launch it by running the file pycalc.pyw. This requires the installation of packages PyQT, numpy and matplotlib. 
If you have the pyinstaller package, you can create an executable version, by running install.bat.

Examples commands:

Numerical expressions: 3+1-4*2    45+6*32*{4[4-[7^2]]}**3    5!*e^(i*pi)    real([4+i]^2), imag[sqrt(i)]    (1,3..9)^2    sum(1..10)    mul([12,10...2])
Fractions: Q(3, 4)+Q(1, 2)   Q[1,Q(3, 4)+Q(1, 2)]    Q(1,6)/Q(1,2)    Q[Q(1,6),Q(1,2)]
Long division: 3x+2x**2+5+4x**4:x+1
Factorization, greatest common divisor and lowest common multiple: decomp(12366)   decomp(x^4 + 6x^3 + 11x^2 + 6x)   parents(24)   parents(24,2)   parents(x^3 - 6x^2 + 11x - 6)   gcd(60,90,450)   lcm(60,90,450)
Algebraic and general expressions: (x+2)^2/(x+2)*(x-4)**2   ### x+1      (### x+1) - {[(x+1)^2]^2}^2    sin(pi*x)exp(-x^2)*e*y^2*sin(pi*x)^3    (sum{(-1)^(k)*x^(2k+1)/(2k+1)! |k=0..n} |n:=0..4), sin(x) from -pi to pi
Algebraic equations and inequalities in one variable: 12x+2=30   (x+3)/5=(x-4)/2   x^5+3x^2=7x-4    x**8=1    x^4 - 10x^3 + 35x^2 - 50x >=-24
Systems of linear equations: 3x+2y-z=5;2x+y+z=6;x+y+z=7    x+y=5;2x+2y=10    2x+y=5;3x+y=7;4x+y=1
Evaluation of general expressions: xy*tan(z) where x=4;z=7     x^2sin(y) |y=pi/2;x=1,3...9     f(x) |f->sin,cos,tan,atan;x=(0,pi/6,pi/4,pi/3,pi/2)     x^2 where x=sin(y) where y=cos(z)     x^2 | x=sin(y) | y=cos(z)     exp(-x)' | x=4    sum([x^2 where x=1,3..11])    mul([x^2 where x=2,4..10])
Ordinary and special functions:sin(x), cos(x), tan(x),asin(x), exp(x), ln(e),log(10),cosh(x),atanh(x) |x=pi/6   sin(pi/[1,2,3,4,6])   sqrt(-1),asin(10),ln(-1)   sqrt(-1+0i),asin(10+0i),ln(-1+0i)    gammafunc(4.5),erf(x),erfc(x),betafunc(4,2), J1(n, x), J2(n, x) | n=0;x=pi/6
Custom functions: f(2,4) where f(x,y):=x*y^2     f(x,g(x)) where f(x,y):=xy and g(x):=sin(x)     sin(x),usin(x) where usin(x):=max(0,sin(x)) from 0 to 8pi   f(x,g(x))|f(x,y):=xy;g(x):=sin(x)    f(x)' |f(x):=sin(x^2),xexp(-x)       q(xy) |q(a):=grad(3a^2),rot(0,0,a)
A few combinatoric and statistical functions:(choose(n,k) |k=0..n) |n=0..5    sum(choose(n,k) |k=0..n) |n=0..8    {binom(p,n,k) where p=0.7;k:=0..n} where n=50   poisson(k,λ) |λ=35;k=0..50    normal(x, μ, σ) |μ=10;σ=1 from 5 to 15    lognormal(x, μ, σ) |μ=1;σ=1 from 0.1 to 10   gammadist(x, α, β) |α=10..30;β=1/2 from 0 to 100   betadist(x, α, β) |α=5;β=1..20 from 0 to 1   MB(v, T, M) |T=273.15;M=29 from 0 to 1500   BT(rate,sensitivity,specificity) |rate=0.15;sensitivity=0.85;specificity=0.9
Ordinary and partial derivatives: sin(cos(x))', (xy^2z^3)_x , (xy^2z^3)_(x,y,z,z)   a^x'    der[sin(x)cos(xy),x,y]
Vector calculus: dot([a, b], [x, y])    [a, b]*[x, y]    [a, b]_1*[x, y]_2    grad(xy)    (y,x^2,z^2x)_x    rot{cross[(y, x, z), (x^2, z^3, y^4)]} |y=1    matmul(R,v) |R={[cos(f),sin(f)],[-sin(f),cos(f)]} where f=pi/4 ;v=(x,y)
Discrete Fourier Transform:FT([1,2,3,4,5])     FT(exp(2πn/20) |n=1..20)    FT(exp(2πin/20) |n=1..20)    FT(1 |n=1..20)    IFT(FT([1,2,3,4,5]))    
1-D ordinary and polar plots: x, x^2, x^3    x', x^2', x^3' data    sin(t) style . pden 100  exp(it)    J1(k, x) |k=1..5 from 0,-1 to 20,1 binom(p, n, k) |n=20, k=12 from 0 to 1 grid on   (x/a)^2+(y/b)**2=1 |a=0.5;b=0.25  sin(4f)  sqrt(cos(2x)) polar pden 20000
1-D parametric plots: (8sin[7t], 3cos[12t])  (sin(t)sqrt{abs[cos(t)]}/[sin(t)+7/5]-2sin(t)+2,t+npi/2) |n=1..4 polar legend off   (sin(2t)^2,cos(9t)+npi/2) |n=0...3 polar   ([b(k-1)]sin(t)-bsin[t(k-1)],[b(k-1)]cos(t)+bcos[t(k-1)]) |b=2;k=0.5,2.5,5 from -6pi to 6pi
2-D plots of functions, vector fields and complex functions: sin(kxy) |k=1,5..20     grad(sin(xy)) from -pi to pi    grad(sin(xy)) from (-pi,-2pi) to (pi,2pi)    (xy,xy^2,xsin(y)^2,cos(xy)) slice 2,4     (r, f)    (x+iy+0.5i)*(x+iy-0.4-0.5i) cont 50 sqrt[(x+iy+0.5)*(x+iy-0.5+0.1i)] cont 0    sqrt[(x+iy+0.5)*(x+iy-0.5+0.1i)] cont 50 reim    3exp(idot(k,x)) |k=(1,4);x=(x,y) reim
Radiation patterns: matmul((-1)^a*r*f+a*f*r,r)/norm(r)^2|a=0,1;f:=(1,0);r:=(x,y)   matmul((-1)^a*r*r/norm(r)^2+a*Id(3),matmul(s*n+n*s,r))/norm(r) |a=0,1; s:=(1,0,0);n:=(0,0,1);r:=(x,y,z);y=0 slice 1,3   matmul{(3r*p-p*r),r}/norm(r)^2 |p:=(0,1);r:=(x,y)
Systems of ordinary differential equations, numerical integration: ode(-y,k) |k=1..10 |k=1..10    ode([dx,dy,0,-g],[0,0,vy,31]) |g=9.8;vy=0,2..10    ode(t^n,0) |n=1..5 t0=0 t1=1
Boundary value problems: bvp([dx,dy,0,-10],[0,0,end,end]) |end=10,20..50   bvp([dx,dy,0,-10],[(1,0,0),(2,0,0),(1,1,xend),(4,1,vend)]) |xend=10;vend=20,10..-20
Iterated maps, bifurcation diagrams: map(-x/1.04,init) |init=1,2,3 tden 120    map(BT(x,0.7,0.5),0.2) data    map([1+u{xcos(t)-ysin(t)},u{xsin(t)+ycos(t)}],[0.9,0.9]) |u=0.9;t=0.4-6/(1+x**2+y**2) tden 2000    map(0.4+0.2i+zexp(kiabs(z)**2),0.2+0.2i) |k=5.1 tden 1000    bifurc[(x,ry(1-y)),(r,0.4)] |r=0,0.005..4 tden=500
Some formulas and simplifications: simp(expr) |expr->sqrt(x)^2,sin(2x),cosh(2x),xsin(x)^2+xcos(x)^2+sinh(x)^2-cosh(x)^2,tan(2x),sinh(2x),exp(ln(x)),simp(exp(ix))


All parentheses are interchangable.

You can use "\" before certain expressions, to convert them to special symbols. For example, writing "\sqrt" and pressing "space" converts "\sqrt" to "√[]". As another example, "\nabla"+"space" and "\dot"+"space" 
yields "∇·", which you can use to calculate the divergence of a vector field. You can also use this to convert Greek letters to their respective symbols (e.g. "\alpha" -> α)

The "|" and "where" symbols are equivalent, they are used to manipulate the expression to their left according to assignments to their right. They can be used inside expressions and nested, 
their scope is limited by the first enclosing pair of parenthesses.
The assignments are performed with the operators "=",":=" and "->". Even though these operators are often interchangable in practice,  their functions differ. 
The "=" operator evaluates the expression at the very end of the calculation, e.g. "x^2' |x=1" first calculates the derivative "2x" of the function "x^2" and then evaluates the result at x=1, yielding 2. 
In contrast, the ":=" operator first assigns the value and then evaluates the expression, so that "x^2' |x:=1" returns 1^2'=0. F
inally, the operator "->" simply replaces all instances of the substituted string with the string to its right, without evaluating it. For example, "fn(x)fnh(x) |f->si" is the same as writing "sin(x)sinh(x)".

There is also a slight difference between the expressions "4,5,6" and "(4,5,6)". The former is to be understood simply as a list of separate expressions, while the latter is more akin to a mathematical vector, 
with the components representing a single entity. For example, asking the parser to plot "t^2,t^3" produces plots of two separate functions of t, while "(t^2,t^3)" produces a parametric curve with the left and 
right expression representing its x and y coordinates, respectively. Similarly, "xy,sin(xy)" produces two 2-D plots while "(xy,sin(xy)" produces a single plot of a 2-D vector field. 
The expression "dot(x,y) |x=(1,2);y=(2,4)" behaves as expected, but with "dot(x,y) |x=1,2;y=2,4" the parser tries to evaluate dot(1,2),dot(1,4),dot(2,2),dot(2,4) and so fails, because the arguments are not vectors.

To plot a function (e.g. sin(x) ), you can either write "plot sin(x)", or you can write "sin(x)" and press ctrl+p. If you enclose the 
To plot a function of two variables (e.g. sin(xy)), you can either write "plot2 sin(xy)", or write "sin(xy)" and press ctrl+i. If the function is complex, (x^2+ixy), the plots of its magnitude and phase will be produced.
If you want to see the real and imaginary part instead, add "reim" to the end of the target expression (e.g. (x^2+ixy reim)).
If the expression is a vector field of two variables, it will be plotted with the arrows, e.g. (x,y^2). If it is a 2-D vector field of one variable (e.g. [cos(t),sin(t)] ), using the "plot" command  (or "ctrl+p") produces a para-
metric plot (circle in this case), while using the "plot2" command (or "ctrl+i") produces a 2-D vector field, assuming that the variable represents the x-coordinate.
If you use variable names such as "r", or "φ", a polar plot will be produced.  You can turn this off by writing "polar off". Conversely, you can convert any plot to a polar plot by writing "polar on".
If you want to animate a function, add anim to the end of the expression to be plotted. For example, you can write "plot sin(tx) anim", or write "sin(tx) anim" and press ctrl+p. 
This also works for functions of three variables, e.g. writing "sin(txy) anim" and pressing ctrl+p produces a 2-D animation. It is always assumed that the  variable that increases with time is the one with lowest rank in the alphabet.
To control the animation speed, use arrows, to pause the animation, press "space" or left-click the expression. Right-clicking the expression reverses the flow of time.
To limit  the time variable from 0 to 10, you can add "t0=0 t1=10". If you want one second of real time to correspond to one second of the mathematical time, add "realtime"
(e.g., the animation "sin(t) t0=0 t1=10 realtime" would last 10 s). 

