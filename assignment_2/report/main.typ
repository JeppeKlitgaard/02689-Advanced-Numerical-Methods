#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "@preview/cheq:0.2.2": checklist
#import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1"
#import "@preview/algorithmic:0.1.0"
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT
// #import algorithmic: algorithm

#import "preamble.typ": mref, mathformatter

#show: checklist
#show: codly-init.with()

#codly(languages: codly-languages)

#set math.equation(numbering: "(1)")

#set heading(numbering: "1.a.i) ")
// #show heading
// #show <nonumber>: set heading(numbering: none)

#set page(
  numbering: "1/1",
  margin: (
    top: 1.5cm,
    bottom: 1.5cm,
    left: 1.5cm,
    right: 2.5cm,
  )
)

#set text(
  size: 11pt
)

#let appendix(body) = {
  set heading(numbering: "A.1", supplement: [Appendix])
  counter(heading).update(0)
  body
}

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: false, upright: true)
#set math.mat(delim:"[")

#let wider = h(3em)

#align(center,
  [#text(
  18pt,
  weight: "bold",
  [Assignment 2])
  #v(-6pt)
  #text(
    14pt,
    [02689 — Advanced Numerical Methods for Differential Equations
  ]
  )
])

#v(8pt)
- Jeppe Klitgaard <`s250250@dtu.dk`>
- Tymoteusz Barcinski <`s221937@dtu.dk`>
- Pernille Christie <`s204249@dtu.dk`>

= Boundary Value Problems <sec:1>

== Legendre Methods <sec:1.a>
We are given the second-order boundary value problem on the domain $x ∈ [a, b]$ where $a=0, b=1$:
$
-ϵ dv(, x, 2)u - dv(, x) u = 1 wider u(a) = u(b) = 0
$ <eq:1_bvp>

Where $ϵ > 0$ is a small parameter characterising the boundary layer width. We are given an exact solution to @eq:1_bvp:
$
u(x) = (e^(-x/ϵ) + (x - 1) - e^(-1\/ε) x)/(e^(-1\/ϵ) - 1)
$ <eq:1_bvp_solution>

=== Confirmation of Exact Solution <sec:1.a_confirm>
We confirm that @eq:1_bvp_solution is a solution to @eq:1_bvp by computing the derivatives:
$
dv(, x)u &= 1/(e^(-1\/ϵ) -1) [ -1/ϵ e^(-x\/ϵ) + 1 - e^(-1\/ϵ)]\
dv(, x, 2)u &= 1/(e^(-1\/ϵ) -1) [ 1/ϵ^2 e^(-x\/ϵ)]\
$

Plugging this into @eq:1_bvp:
$
-ϵ dv(, x, 2)u - dv(, x) u = 1/(e^(-1\/ϵ) -1) [cancel(-ϵ 1/ϵ^2 e^(-x\/ϵ) + 1/ϵ e^(-x\/ϵ)) - 1 + e^(-1\/ϵ)] = 1
$
We observe that the proposed solution indeed satisfies @eq:1_bvp.

=== Coordinate Transformation

The implementation of of numerical schemes such as the Spectral Legendre Tau Method (LTM) and the Legendre Collocation Method (LCM) are typically undertaken in the canonical domain $z ∈ [-1, 1]$, which differs from our problem domain $x ∈ [0, 1]$. This can be remedied by either performing a coordination transformation of _differential operator_ of the problem into the canonical domain, or alternatively by transforming the _differentiation matrix_ of the method.

In either case, we consider the coordinate transformation at hand, with associated differentials:
$
x = (b - a)/2 z + (a+b)/2, quad dif x = (b - a)/2 dif z
$ <eq:1_coordination_transformation>

From which we can compute the change of basis using the chain rule:
$
dv(,x) &= dv(z, x) dv(,z) = 2/(b-a) dv(,z)\
dv(,x,2) &= dv(z, x, 2) dv(,z) + (dv(z, x))^2 dv(,z,2) = 4/(b-a)^2 dv(,z,2)
$

Opting to transform our problem into the canonical domain results in the following modification of @eq:1_bvp for problem domain $x ∈ [a, b] = [0, 1]$:
$
-4ϵ dv(, z, 2) u(z) - 2 dv(, z) u(z) = 1 wider u(-1) = u(1) = 0
$ <eq:1_bvp_canonical>

With associated analytical solution based on @eq:1_bvp_solution:
$
u(z) = (e^(-(z + 1)/(2ϵ)) + (z-1)/2 - e^(-1\/ε) (z+1)/2)/(e^(-1\/ϵ) - 1)
$ <eq:1_bvp_solution_canonical>

=== Spectral Legendre Tau Method (LTM) <sec:1.a_ltm>

The _Spectral Legendre Tau Method_ may be formulated as a _Method of Weighted Residuals_ (MWR)
wherein the residual of the $N$th order approximation of a function $u(z)$ is given by:
$
cal(R)_N (z) = u(z) - u_N (z) = cal(L) u_N (z) - f(z)
$

Where $cal(L)$ and $f(z)$ are the _differential operator_ and _right-hand side function_ respectively, while $u_N (z)$ is the $N$th order approximation as given by:
$
u_N (z) = sum_(n=0)^N hat(u)_n ϕ_n (z)
$ <eq:1_approx_expansion>

With ${ϕ_n (x)}_(n=0)^N$ being orthogonal basis functions defined by their inner product and weight $w(z)$:
$
(ϕ_i, ϕ_j)_w = γ_i δ_(i j)
$

The MWR family of methods are characterised by the projection of the residual $cal(R)_N$ onto a set of _test functions_ ${ψ_n}_(n=0)^N$ vanishing, where the the choice of test functions determine whether the method is the Collocation Method, Tau Method, or some other method.

For the Tau Method, the test functions are given by the basis functions themselves, such that $ψ_n (z) = ϕ_n (z)$, which by the MWR condition leads to:
$
(cal(R)_N, ϕ_i)_w = ∫_a^b (u - sum_(n=0)^N hat(u)_n ϕ_n) ϕ_i w dif x = 0 wider i = 0, 1, ..., N
$

And orthogonality leads us to an expression for the expansion coefficients $hat(u)_n$:
$
hat(u)_n = (u, ϕ_n)_w / (ϕ_n, ϕ_n)_w wider n = 0, 1, ..., N
$

However, unlike the _Galerkin Method_, the Tau Method relaxes the requirement that the basis functions ${ϕ_n (x)}_(n=0)^N$ satisfy the boundary conditions. As such, we must _spend_ two of the $N+1$ independent equations we were granted by @eq:1_approx_expansion to abide by the boundary conditions, leaving us with $N-1$ remaining equations to constrain using the projection of the residual onto the basis functions:
$
sum_(n=0)^N hat(u)_n (cal(L) ϕ_n, ϕ_i)_w = (f, ϕ_i)_w = hat(f)_i wider i = 0, 1, ..., N-2
$
Where the number of indices $i$ is determined by the remaining number of equations.

The boundary conditions are imposed using the remaining two equations:
$
cal(B)_- u_N (α) = g_- wider cal(B)_+ u_N (β) = g_+ wider
$
Where $cal(B)_±$ are the differential operators associated with the conditions (most often _Dirichlet_ or _Neumann_) while $α, β$ are the boundary points and $g_±$ are the boundary values.

// JK: Above is a general introduction with symbols that we can easily refer to below

From this, we may construct a linear system of equations, $cal(L)_N$:
$
cal(L)_N hat(vv(u)) = hat(vv(f))
$ <eq:1_lin_sys_1>

Which gives modal solution:
$
hat(vv(u)) &= (cal(L)_N)^(-1) hat(vv(f))\
u_N (x) &= sum_(n=0)^N hat(u)_n ϕ_n (x)
$

In order to find the discretisation given by @eq:1_lin_sys_1, we may obtain a sparse matrix $mm(A)$ by finding a recurrence relation dictating the coefficients of the algebraic equations.

For this exercise, we use the Legendre basis functions, i.e.
$
u(z)=sum_(n=0)^infinity hat(u)_n L_n (z) quad quad dv(,z,q)u(z)=sum_(n=0)^infinity hat(u)_n^((q))L_n(z), 
$
where the coefficients for the derivatives are given by the recurrence
$
hat(u)_n^((q-1))= hat(u)_(n-1)^((q))/(2n-1)-hat(u)_(n+1)^((q))/(2n+3), quad n>=1.
$ <eq:1_recursion>
Furthermore, the derivative coefficients are related to the solution coefficients through
$
hat(u)_n^((1))=(2n+1)sum_(p=n+1\ n+p "odd")^infinity hat(u)_p, wider hat(u)_n^((2))=(n+1/2)sum_(p=n+2\ n+p "even")^infinity (p(p+1)-n(n+1))hat(u)_p, quad n>=0.
$<eq:1_uhat_expression>
Using the recurrence, we can build the system matrix $A$.

// note that we need to apply coordinate transformation but it cancels out on both side of the weak form

//weight function for Legendre w(x) = 1, so convenient

Consider a general second-order equation and its rewritten version,
$
  a hat(u)^((2))_n + b hat(u)^((1))_n + c hat(u)_n= hat(f)_n, quad hat(u)^((2))_n =
  -(b)/(a) hat(u)^((1))_n - c/a hat(u)_n + 1/a hat(f)_n.
$<eq:1_eq_sys>
Using the recursion for derivatives of Legendre-polynomials, @eq:1_recursion, we get the following for $q=2$ and $q=1$
$
  hat(u)^((1))_n = 1/(2n-1) hat(u)^((2))_(n-1) - 1/(2n+3) hat(u)^((2))_(n+1), quad  hat(u)_n = 1/(2n-1) hat(u)^((1))_(n-1) - 1/(2n+3) hat(u)^((1))_(n+1).
$
We plug this into @eq:1_eq_sys
$
  hat(u)^((1))_n = 1/(2n-1) (-(b)/(a) hat(u)^((1))_(n-1) - c/a hat(u)_(n-1) + 1/a hat(f)_(n-1)) - 1/(2n+3) (-(b)/(a) hat(u)^((1))_(n+1) - c/a hat(u)_(n+1) + 1/a hat(f)_(n+1)),
$
and rearrange the terms,
$
    hat(u)^((1))_n = &-(b)/(a) (1/(2n-1) hat(u)^((1))_(n-1) - 1/(2n+3) hat(u)^((1))_(n+1) ) - c/a (1/(2n-1) hat(u)_(n-1) - 1/(2n+3) hat(u)_(n+1) ) \
    &+ 1/a (1/(2n-1) hat(f)_(n-1) - 1/(2n+3) hat(f)_(n+1))
$
We can now use the recursion for $q=1$ to simplify,
$
    hat(u)^((1))_n = -(b)/(a) hat(u)_n - c/a (1/(2n-1) hat(u)_(n-1) - 1/(2n+3) hat(u)_(n+1) ) + 1/a (1/(2n-1) hat(f)_(n-1) - 1/(2n+3) hat(f)_(n+1)).
$
Consider now
$
  1/(2n-1) hat(u)^((1))_(n-1) =& 1/(2n-1) [-(b)/(a) hat(u)_(n-1) - c/a (1/(2n-3) hat(u)_(n-2) - 1/(2n+1) hat(u)_(n) ) + 1/a (1/(2n-3) hat(f)_(n-2) - 1/(2n+1) hat(f)_(n))]  \
  =& (- c/a 1/(2n-1) 1/(2n-3)) hat(u)_(n-2) + (- b/a 1/(2n-1)) hat(u)_(n-1) + (c/a 1/(2n-1) 1/(2n+1)) hat(u)_(n) \
  &+ (1/a 1/(2n-1) 1/(2n-3)) hat(f)_(n-2) + (- 1/a 1/(2n-1) 1/(2n+1)) hat(f)_(n),

$<eq:1_uhat1_term1>
and likewise
$
  1/(2n+3) hat(u)^((1))_(n+1) =& 1/(2n+3) (-(b)/(a) hat(u)_(n+1) - c/a (1/(2n+1) hat(u)_(n) - 1/(2n+5) hat(u)_(n+2) ) + 1/a (1/(2n+1) hat(f)_(n) - 1/(2n+5) hat(f)_(n+2)))  \
  =& (- c/a 1/(2n+3) 1/(2n+1)) hat(u)_(n) + (- b/a 1/(2n+3)) hat(u)_(n+1) + (c/a 1/(2n+3) 1/(2n+5)) hat(u)_(n+2) \
  &+ (1/a 1/(2n+3) 1/(2n+1)) hat(f)_(n) + (- 1/a 1/(2n+3) 1/(2n+5)) hat(f)_(n+2).
$<eq:1_uhat1_term2>

We can now use the recursion for $q=1$ once again,
$
  1/(2n-1) hat(u)^((1))_(n-1) - 1/(2n+3) hat(u)^((1))_(n+1) - hat(u)_(n) = 0
$ <eq:1_recursion_rewritten>
Using @eq:1_uhat1_term1 and @eq:1_uhat1_term2 in @eq:1_recursion_rewritten, we get an equation on the form
$
  a_(n, n-2) hat(u)_(n-2) + a_(n, n-1) hat(u)_(n-1) + a_(n, n) hat(u)_(n) + a_(n, n+1) &hat(u)_(n+1) + a_(n, n+2) hat(u)_(n+2) = g_(n-2) hat(f)_(n-2) + g_(n) hat(f)_(n) + g_(n+2) hat(f)_(n+2). quad
$
Clearly, there are overlapping terms for $hat(u)_(n)$ and $hat(f)_(n)$. For the $hat(u)_n$ terms, we get
$
  a_(n,n)hat(u)_n=&(c/a 1/(2n-1) 1/(2n+1)) hat(u)_(n) - (- c/a 1/(2n+3) 1/(2n+1)) hat(u)_(n) - hat(u)_(n) \
  =& (c/a 1/(2n-1) 1/(2n+1) + c/a 1/(2n+3)1/(2n+1) - 1) hat(u)_(n) \
  =& ( c/a 2 / ((2n-1)(2n+3)) - 1)hat(u)_(n),
$
and for the $hat(f)_n$ terms,
$
  g_n hat(f)_n &=(- 1/a 1/(2n-1) 1/(2n+1)) hat(f)_(n) - (1/a 1/(2n+3) 1/(2n+1)) hat(f)_(n) \
  &= - (1/a 1/(2n-1) 1/(2n+1) + 1/a 1/(2n+3) 1/(2n+1)) hat(f)_(n) \
  &=
  - 1/a 2 / ((2n-1)(2n+3)) hat(f)_(n).
$
Thus, the coefficients for the linear system are the following
$
  a_(n, n-2) &= - c/a 1/(2n-1) 1/(2n-3), \
  a_(n, n-1) &= - b/a 1/(2n-1), \
  a_(n, n)   &= c/a 2 / ((2n-1)(2n+3)) - 1, \
  a_(n, n+1) &= b/a 1/(2n+3),\
  a_(n, n+2) &= - c/a 1/(2n+3) 1/(2n+5),\
  g_(n-2) &= - 1/a 1/(2n-1) 1/(2n-3),\
  g_(n) &= 1/a 2 / ((2n-1)(2n+3)),\
  g_(n+2) &= - 1/a 1/(2n+3) 1/(2n+5).\
$<eq:1_tau_system>
These coefficients are of course only valid for $2 <= n < N-2$, and we need $N$ equations to solve the system uniquely. To get two more equations, we utilize @eq:1_uhat_expression. We are of course not able to calculate the infinite series, so we use a truncated version for $n=0,1$,
$
1/2 a sum_(p=2\ p "even")^(N-1) p(p+1)hat(u)_p+b sum_(p=1\ p "odd")^(N-1) hat(u)_p+c hat(u)_0&=hat(f)_0,\
3/2 a sum_(p=3\ p+1 "even")^(N-1) (p(p+1)-2) hat(u)_p+3b sum_(p=2\ p+1 "odd")^(N-1) hat(u)_p + c hat(u)_1&=hat(f)_1.
$
For the remaining two needed equations, we impose the boundary conditions, namely
$
sum_(n=0)^(N-1) hat(u)_n L_n (-1)&=0,\
sum_(n=0)^(N-1) hat(u)_n L_n (1)&=0.
$

With these equations, we set up a new system, $cal(vv(A))hat(vv(u))=hat(vv(g))$, where the coefficients are given by @eq:1_tau_system for $2 <= n < N-2$ and 
$
a_(0,n)& = cases(
  c &", " n=0,
  1/2 a n(n+1) &", " n>0 "and even" ,
  b &", " n>0 "and odd"  
)\
a_(1,n) &= cases(
  0 &", " n=0,
  c & ", " n=1, 
  3/2 a (n(n+1)-2) &", " n>1 "and odd" ,
  b &", " n>1 "and even"  
)\
a_(N-2,n)&=L_n (-1)\
a_(N-1,n)&=L_n (1)\
hat(g)_0&=hat(f)_0\
hat(g)_1&=hat(f)_1\
hat(g)_(N-2)&=0    \
hat(g)_(N-1)&=0    
$


=== Legendre Collocation Method (LCM) <sec:1.a_lcm>

The Legendre Collocation Method is somewhat simpler to implement and derive.

As it is a collocation method, we employ the _Dirac Delta_ function as our test functions rather than the basis functions themselves as was done for the Tau Method:
$
ψ_i (x) = δ(x - x_i)
$

By the MWR condition, we require the residual $cal(R)_N$ to disappear at every nodal point $x_i$:
$
(cal(R)_N, δ(x - x_i))_w
&= ∫_a^b cal(R)_N (x) δ(x-x_i) w(x) dif x = 0\
&= cal(R)_N (x_i) = sum_(n=0)^N hat(u)_n ϕ_n (x_i) - u(x_i)
$
Where we have used $w(x) = (1-x)^α (1+x)^β = 1$ for the Legendre polynomials characterised by $α=β=0$. This yields the collocation condition:
$
u(x_i) = sum_(n=0)^N hat(u)_n ϕ_n (x_i) wider i = 0, 1, ..., N
$ <eq:1_collocation>

Where we may determine the expansion coefficients ${hat(u)_n}_(n=0)^N$ by orthogonality:
$
hat(u)_N = ((u, ϕ_n)_(w, N))/((ϕ_n, ϕ_n)_(w, N)) = 1/γ_n ∫_a^b u_N (x) ϕ_n (x) dif x wider n = 0, 1, ..., N
$ <eq:1_modal_expansion_coeffs>

Where the integral will usually be computed using a quadrature. If we are free to pick our collocation points, $x_i$, and quadrature rule, we may pick the Legendre-Gauss-Lobatto points and a Jacobi Gauss-type quadrature, which is exact up to a polynomial degree of $2N+1$ @kopriva[p. 33].

Recalling the construction of the _Vandermonde matrix_, $cal(V)$ from the first assignment, we find that it corresponds exactly to the conversion between the nodal and modal bases outlined above, thus allowing the expression of @eq:1_modal_expansion_coeffs as:
$
u_N = cal(V) hat(u)_N
$

Where:
$
cal(V) = mat(
  ϕ_1 (x_1), ϕ_2 (x_1), ..., ϕ_N (x_1);
  ϕ_1 (x_2), ϕ_2 (x_2), ..., ϕ_N (x_2);
  ⋮, ⋮, ⋱, ⋮;
  ϕ_1 (x_N), ϕ_2 (x_N), ..., ϕ_N (x_N);
)
$

This in turn enables the construction of the _differentiation matrix_ with respect to the canonical coordinate $z$, $cal(D)_z$:
$
cal(D)_z = cal(V)_z cal(V)^(-1)
$

Where $cal(V)_z$ is constructed in similar fashion to $cal(V)$, though using the derivatives of the basis functions with respect to $z$.

Using the differentiation matrix $cal(D)_z$ we are able to construct a discrete operator for the differential equation given in @eq:1_bvp_canonical:
$
  cal(L)_N u(z) = -4ϵ cal(D)_z^2 u(z) - 2 cal(D)_z u(z) = f(z) = 1
$

Lastly, the boundary conditions are imposed by replacing the first and last rows of the operator $cal(L)_N$, corresponding to the boundary points $z = -1$ and $z = 1$, with a mask that selects the corresponding element of $u_N$ while updating the right-hand side to the boundary values accordingly.

=== Results <sec:1_a_resluts>
To check out solvers, we first try to solve it and compare it to the exact solution. The plots for the Legendre Tau Method can be seen in @fig:1_sol_tau_eps0001, @fig:1_sol_tau_eps001 and @fig:1_sol_tau_eps01. 
#figure(
  image("output/1_a_tau_eps0001.png"),
  caption: [Computed solutions using LTM for different values of $N$ to the problem with $epsilon=0.001$]
)<fig:1_sol_tau_eps0001>
#figure(
  image("output/1_a_tau_eps001.png"),
  caption: [Computed solutions using LTM for different values of $N$ to the problem with $epsilon=0.01$]
)<fig:1_sol_tau_eps001>
#figure(
  image("output/1_a_tau_eps01.png"),
  caption: [Computed solutions using LTM for different values of $N$ to the problem with $epsilon=0.1$]
)<fig:1_sol_tau_eps01>

Likewise, he plots for the Legendre Collocation Method can be seen in @fig:1_sol_col_eps0001, @fig:1_sol_col_eps001 and @fig:1_sol_col_eps01. 
#figure(
  image("output/1a_plot_working.png"),
  caption: [Computed solutions using LCM for different values of $N$ to the problem with $epsilon=0.001$]
)<fig:1_sol_col_eps0001>
#figure(
  image("output/1_a_col_eps001.png"),
  caption: [Computed solutions using LCM for different values of $N$ to the problem with $epsilon=0.01$]
)<fig:1_sol_col_eps001>
#figure(
  image("output/1_a_col_eps01.png"),
  caption: [Computed solutions using LCM for different values of $N$ to the problem with $epsilon=0.1$]
)<fig:1_sol_col_eps01>

For both methods, we see that we need quite different grid refinements to achieve results that mimic the solution. To further investigate this, we consider the errors for each of the methods.


=== Discussion of Methods <sec:1.a_discussion>

We may evaluate the performance of our models above by considering the relative error against the known solution given by @eq:1_bvp_solution by again using a quadrature to compute a highly accurate estimate of the $L^2$ norm over the domain $[0, 1]$:
$
norm(u_N (x) - u(x))_(L^2_[0, 1]) ≈ sqrt(sum_(j=0)^N (u_N (x_j) - u(x_j))^2 J w_j)
$

Where $J = (b-a)/2=1/2$ is the Jacobian from the coordinate transformation given in @eq:1_coordination_transformation and $w_j$ are the weights from the Jacobi-Gauss quadrature that was used to determine the abscissas $x_j$ as detailed in @kopriva[eq. 1.131].

This allows us to investigate the convergence of the methods by plotting the relative error against $N$ for varying values of $epsilon$, see @fig:1_error_collocation and @fig:1_error_tau.
#figure(
  image("output/1_a_collocation_convergence.png"),
  caption: [Plot of the relative error for a solution using LCM for varying $N$.]
)<fig:1_error_collocation>

#figure(
  image("output/1_a_tau_convergence.png"),
  caption: [Plot of the relative error for a solution using LTM for varying $N$.]
)<fig:1_error_tau>

The plots show that the convergence is quite similar for the two methods and that they perform quite equally. As $epsilon$ grows, the convergence is much faster, getting to machine-precision around $N=25$ for $epsilon=0.01$ compared to $N=175$ for $epsilon=0.001$. To explain this difference, we look at the exact solutions, as plotted, see Figures 1-6. Clearly the solution is much steeper close to $x=0$ for smaller values of $epsilon$. This means that the solver will need a much finer grid to achieve enough precision in the estimation of the derivative close to $x=0$, explaining the much slower convergence rate.

A disadvantage of both methods is that they have even spacing, which means that there is no other way to accommodate for steeper gradients than simply increasing the grid refinement $N$. Consequently, we get the slow convergence rate for small $epsilon$ that we just saw.

The Legendre Tau Method results in a sparse and banded system matrix (excluding the boundary conditions and first 2 rows from @eq:1_uhat_expression ) which means that it is more memory efficient and the system of equations can be solved in linear time with respect to $N$.  This is of course a large advantage, when we consider problems, such as the problem presented here for $epsilon=0.001$, where a large number of grid points is needed to get precise solutions. Note that the sparsity of the system matrix is obtained by using the recursions and implementing the system matrix only with @eq:1_uhat_expression results in dense upper triangular matrix.

== Irrotational Flow Around a Cylinder <sec:1.b>

We are given the Laplace's Equation applied to the scalar velocity potential $ϕ(x, y)$ in 2 dimensions expressed in polar coordinates:
$
∇^2 ϕ = 1/r pdv(, r) (r pdv(ϕ, r)) + 1/r^2 pdv(ϕ, θ, 2) = 0, wider (r, θ) ∈ [r_1, ∞[ × [0, 2π], wide r_1 > 0 
$ <eq:1b_diff_eq>

With the exact solution:
$
ϕ(r, θ) = V_∞ (r + r_1^2/r) cos(θ)
$ <eq:1b_exact_sol>

Where $V_∞ = 1$ is the far-field flow speed and $r_1$ is the radius of the obstructing cylinder.

In order to carry out numerical analysis, we require the domain to be finite, to which end we introduce $r_2 ∈ ℝ$ such that $0 < r_1 < r_2$.

For our discretisation, we choose to keep the original $θ$ coordinate, which we solve using the spectral Fourier collocation method in order to exploit the intrinsic periodicity of the angular coordinate.

For the $r$ coordinate, we employ the spectral Legendre collocation method which has canonical domain $z_r ∈ [-1, 1]$, which we may transform into our physical domain using the transformation with associated Jacobian:
$
r(z_r) = (r_2 - r_1)/2 z_r + (r_1 + r_2)/2, quad J_r = (r_2 - r_1) / 2
$
Thus allowing us to re-express the original differential equation, @eq:1b_diff_eq, restricted to a transformed finite radial domain $[-1, 1]$ as:
$
∇^2 ϕ
&= 1/r(z_r) J_r^(-1) pdv(, z_r) (r(z_r) J_r^(-1) pdv(ϕ, z_r)) + 1/r(z_r)^2 pdv(ϕ, θ, 2)\
&= (J_r^(-1))^2 pdv(ϕ, z_r, 2) + 1/r(z_r) J_r^(-1) pdv(ϕ, z_r) + 1/r(z_r)^2 pdv(ϕ, θ, 2)
= 0
$

Where the domain is given by $(z_r, θ) ∈ [-1, 1] × [0, 2π]$.

Construction of the differentiation matrices is similar to the treatment in @sec:1.a_lcm, though with the use of the discrete Fourier series as the basis functions, $ψ(x)$:
$
dv(ψ(x), x, k) = ((i 2π n) / L)^k e^((i 2π n )/L x)
$

With this, we are able to implement the differentiation matrices for the $z_r$ and $θ$ coordinates and use the standard technique of applying the Kronecker product to produce a single system of equations that capture both coordinates:
$
cal(D)_(z_r,"full") &= cal(D)_(z_r) times.o I_θ\ 
cal(D)_(θ,"full") &= I_(z_r) times.o cal(D)_θ\ 
$

Where $I_(|⋅|)$ are the identity matrices associated with the respective coordinates. 

With this, we are able to construct the following discretisation of the problem, implicitly letting $r = r(z_r)$ and using the '$"full"$' subscript to denote the construction that matches the multi-dimensional differentiation operators:
$
cal(L)_N ϕ = (J^(-1))^2 cal(D)_(z_r,"full")^2 ϕ + 1/r_"full" J^(-1) cal(D)_(z_r,"full") ϕ + 1/r^2 cal(D)_(θ,"full")^2 ϕ
$

After which the solution $ϕ$ may be computed as:
$
ϕ_N = (tilde(cal(L))_N)^(-1) tilde(f)_N
$

Where the overhead tildes denote that the incorporation of boundary conditions on the obstructing cylinder and the domain boundary as given by $ϕ_"cylinder"(θ) = ϕ(r_1, θ)$ and $ϕ_"boundary"(θ) = ϕ(r_2, θ)$ respectively.

Solving this yields a solution for the scalar velocity potential $ϕ(r, θ)$, which in turn can be used to compute a velocity field $(u, v) = ∇ϕ$, which may be represented alongside the velocity potential as seen below in @fig:1b_main, in which we observe the expected motion of the fluid around obstructing cylinder.

#figure(
  image("output/1_b.png"),
  caption: [
    Solution to irrotational flow around an obstructing cylinder with both exact and numerical solutions shown. 
  ]
) <fig:1b_main>

We then investigate the convergence behaviour of our implementation as shown in @fig:1b_convergence. Notably, we find that the solver is _exact_ in the angular direction for angular discretisation $N_θ > 1$, which may be understood by referring to the exact solution from @eq:1b_exact_sol:
$
ϕ(r, θ) = V_∞ (r + r_1^2/r) cos(θ)
$

Clearly the angular dynamics will be captured exactly as soon as the first order trigonometric polynomial, $e^(i θ)$ is included, which exactly captures the $cos(θ)$ factor that constitutes the entire angular behaviour of the solution.

We may understand the behaviour for the $N_θ=1$ case by asserting that the error is dominated by the angular discretisation in this case. Importantly, $N_θ$ denotes the _discretisation order_, which is one higher than the polynomial order, $P_θ=N_θ-1$. As such, the $N_θ=1$ case corresponds to a single mode capturing the signal mean only.

For $N_θ>1$ the error is dominated by the radial discretisation and we observe the expected spectral convergence, which achieves machine precision at around $N_r = 20$.

#figure(
  image("output/1_b_convergence.png"),
  caption: [
    Convergence plot for the irrotational flow around a cylinder solved using Legendre polynomials in the radial direction and trigonometric polynomials in the angular direction.
  ]
) <fig:1b_convergence>


= Time-Dependent Problems <sec:2>

#counter(heading).update((2, 2))



== Spectral Fourier Solver for the KdV Equation <sec:2.c>
We are given the following non-linear initial value problem
$
  partial / (partial t) u + 6 u partial / (partial x) u + partial^3 / (partial x^3) u = 0, quad -infinity < x < infinity, quad t > 0
$
with an exact solution given as
$
  u(x,t)=f(x-c t), quad f(x)=1/2 c sech^2(1/2sqrt(c)(x-x_0))
$
where $x_0$ is the center of the soliton at $t=0$ and $c>0$ is the speed of the soliton as it travels to the right. With this knowledge of a solution, we can calculate an initial condition, to be used in our solver,
$
  u(x, 0) = f(x).
$
To solve this problem, we will treat space and time separately using Method of Lines. For the spatial discretisation, we use the spectral Fourier Collocation method. Along the spectral dimension $x$ we have 
$
  lim_(x -> infinity)u(x,t)=0,
$<eq:2_lim_condition>
for finite $t$. Thus, if we go out far enough, we may approximate the function by the periodic Fourier basis, since the left and right boundaries are approximately equal for all times (since they will be approximately 0). We will therefore consider $x$ on a finite interval $[-alpha, alpha]$ for some large $alpha$. Since the Fourier basis functions are defined on $[0,2pi]$, we have to be careful when evaluating. However, if we simply scale the frequencies by $(2pi)/L$, where $L$ is the length of the interval, when evaluating the Fourier functions, everything will work.

/*
We define the coordinate transformation to map the problem from $[-alpha, alpha]$, where $alpha$ is some big number, to the $[0, 2 pi]$ interval such that the condition in @eq:2_lim_condition is satisfied.
$
  x = ((z - pi)/pi)alpha = alpha / pi z - alpha, quad partial / (partial x) = pi / alpha partial / (partial z).
$
We then rewrite the differential equation,
$
  partial / (partial t) u + (6pi)/alpha u partial / (partial z) u + pi^3 / alpha^3 partial^3 / (partial z^3) u = 0, quad z in [0, 2pi], quad t > 0.
$
and the exact solution exact solution
$
  f(x) = 1/2 c sech^2(1/2 sqrt(c)(x - x_0)), quad u(x, t) = f(x - c t)
$
becomes
$
  tilde(f)(z) = f(alpha / pi z - alpha) = f(x)
$
*/

// COMMENT HERE THAT IT IS EFFECTLIVLY THE COORDIANTE TRANSFORMATION.

Using the aforementioned Fourier basis with scaled frequencies, we can write an expression for the approximated solution
$
  u_N (x,t) = sum_(k=-N\/2)^(N \/ 2) hat(u)_k (t) phi.alt_k (x), quad phi.alt_k (x) = e^(i (2 pi k)/(2alpha) x )=e^(i (pi k)/(alpha) x )
$

We use the Fourier differentiation matrix $D$ constructed from the Vandermonde matrices with the Fourier basis functions.
// from assignment 1 with explicit formulas for Lagrange polynomials, $D$.
We use the spectral Fourier Collocation Method, and we therefore want the residual to be zero at each collocation node,
$
  lr({(u_N)_t+6u_N (u_N)_x + (u_N)_(x x x)}  |)_(x_j) = 0, wider j=0,1,...,N
$
We look at the derivatives of $u_N$,
$
  partial_x u_N (x_j, t)=sum_(k=-N \/ 2)^(N\/2) hat(u)_k (t) i (pi k)/alpha phi_k (x_j), quad partial_(x x x) u_N (x_j, t)=-sum_(k=-N/2)^(N/2) hat(u)_k (t) i((pi k)/alpha)^3 phi_k (x_j),
$
and insert this into the equation. Defining $vv(u)_N (t)=[hat(u)_0 (t), hat(u)_1 (t),..., hat(u)_N (t)]$, we get a nonlinear operator
$
  L_N bold(u)_N (t) = -6 "diag"(vv(u)_N (t)) D vv(u)_N (t) - D^3 vv(u)_N (t)
$<eq:L_N_nonlienar>
We obtain a system of ordinary differential equations with initial condition
$
  partial / (partial t) vv(u)_N (t) = L_N bold(u)_N (t), quad vv(u)_N (0) = [f(x_0), f(x_1),..., f(x_N)].
$
To solve the system of ordinary differential equations we use `solve_ivp` form SciPy's `integrate` module with it's standard numerical solver `RK45`, which is an explicit Runge-Kutta method of order 5(4), @solve_ivp.

To ensure that we get a stable solver, we have to control the time stepping. Specifically, we have to control the size of each time step, $Delta t$. For a linear system of equations of the form
$
  (d y) / (d t) = A y, quad y in R^m, quad A in R^(m times m)
$
the stability condition is the following @L6_slides[slide 40]
$
  Delta t <= s /(C N^p), quad lambda_"max" (A)<=C N^p.
$<eq:stability_region>
where the parameter $s$ is derived from the the absolute stability region of the ODE solver. We have found that for Runge-Kutta 4 $s=2.7853$ @s_RK45.

Since our system matrix is nonlinear, see @eq:L_N_nonlienar, we consider its linearization such that we could apply the above result. To this end, we apply the technique of frozen coefficients from @L7_slides[slide 36], where we approximate the system matrix as
$
  dash(L)_N (t = 0) approx - 6 (max_(i in {0, ... ,N}) |u_N (x_i, t = 0)|) D - D^3,
$
which results in a linear operator $dash(L)_N (0)$ and allows us to find the its eigenvalues and relate them to the absolute stability region of the Runge-Kutta solver and the condition @eq:stability_region.

We note that for the considered KdV problem we have that the $max_(x in [0, 2pi]) |u(x, t)| = c/2$ is constant over time and therefore
$
max_(i in {0, ... ,N}) |u_N (x_i, 0)| approx max_(i in {0, ... ,N}) |u_N (x_i, t)| quad => dash(L)_N (0) approx dash(L)_N (t) quad forall t
$
where the approximation is a result of the discretization. This implies that the eigenvalues $lambda(dash(L)_N (t))$ should not change with time. If it was not the case, then the eigenvalues of the linear operator $lambda(dash(L)_N (t))$ should be computed during the time integration of the ODE solver and the stability condition should be updated accordingly.

Consider again @eq:stability_region. For the $p$ value, we have $|lambda_"max"|=max|i|=N$ for $D$, @L7_slides[slide 37]. In this case, we have $D^3$ in our system, so we must have $p=3$. To find a value for $C$ and confirm the value for $p$, we plot $max_(i in {0, ... ,N}) |lambda(dash(L)_N (0))|$, see @fig:2_max_eigvals.

#figure(
  image("output/2_c_max_eigvals.png"),
  caption: [Plot of $max_N|lambda|$ for varying values of $N$. The equation for the reference line is $C N^p$ with $C=3 dot 10^(-4)$ and $p=3$.]
) <fig:2_max_eigvals>

The plot confirms that $p=3$ and furthermore, we see that choosing $C=3 dot 10^(-4)$ ensures  $lambda_"max" (dash(L)_N (0))<=C N^p$  for all large values of $N$ as we want. 

We can now check the performance of our solver by plotting the computed solutions and compare them to the exact solution for varying $t$. This can be seen in @fig:2_c_solution.

#figure(
  image("output/2_c_solutions.png"),
  caption: [Plot of the numerical solutions and the exact solutions for varying $t$.]
) <fig:2_c_solution>

The solver seems to work as wanted. However, we need to do more analysis of the error, and we will do this in the following questions.


== Testing the Solver <sec:2.d>
To further test the solver, we calculate the estimated errors $||u- cal(I)_N u||_2$ and $||u- cal(I)_N u||_infinity$ for varying values of the soliton speed, $c$, and number of nodes, $N$. We start by doing a plot for $alpha=25$, see @fig:2_d_alpha25.
#figure(
  image("output/2_d_error_alpha25.png"),
  caption: [The estimated $L_2$- and $L_infinity$-errors with $alpha=25$ and for different values of soliton speed, $c$, and varying number of nodes, $N$.]
) <fig:2_d_alpha25>
None of the solutions seem to converge, on the contrary the error grows after $N=50$. This might stem from the fact that we consider the solution in a too small finite interval, which means that our periodic boundary conditions are not satisfied for the exact solution, we are seeking. We therefore try to solve the problem with $alpha=50$, see @fig:2_d_alpha50.

#figure(
  image("output/2_d_error_alpha50.png"),
  caption: [The estimated $L_2$- and $L_infinity$-errors with $alpha=50$ and for different values of soliton speed, $c$, and varying number of nodes, $N$.]
) <fig:2_d_alpha50>
In this case, we achieve convergence for all values of $c$. While we reach machine precision for $c=0.5$ and $c=1.0$, we are only able to reach approximately $10^(-11)$ for both of the estimated errors for $c=0.25$.

The higher floor for a lower $c$ must stem from the fact that the soliton will be more tail-heavy for lower $c$. Thus we need a larger interval to ensure that the periodic boundary conditions are upheld. This means that it might be worth to decide $alpha$ depending on $c$ to ensure that we will always be able to find a numerical solution.


Another way to validate the solver is to check that the three fundamental quantities mass ($M$), momentum ($V$) and energy ($E$) are all invariant with respect to time.
To approximate the mass, momentum and energy, we only consider the interval $[-alpha, alpha]$ to allow us to estimate the integral using the trapezoidal rule as follows
$
M =  integral_(-infinity)^infinity u(x) dif x approx integral_(-alpha)^alpha u(x) dif x approx sum_(i=1)^N (u(x_(i-1))+u(x_i))/2 Delta x_i =  tilde(M),
$
and likewise for momentum and energy,
$
V&=  integral_(-infinity)^infinity u^2(x) dif x approx sum_(i=1)^N ((u(x_(i-1)))^2+(u(x_i))^2)/2 Delta x_i=tilde(V), \
E&=integral_(-infinity)^infinity (1/2 u_x^2(x) +u^3(x)) dif x \
&approx sum_(i=1)^N ((1/2 (u_x (x_(i-1)))^2-(u(x_(i-1)))^3)+((1/2 (u_x (x_(i)))^2-(u(x_i))^3)))/2 Delta x_i=tilde(E).
$
In the code, we use the function `trapezoidal` from SciPy's `integrate` module. The approximated quantities over time can be seen in @fig:2_d_conservation for varying $N$.

#figure(
  image("output/2_d_conservation.png"),
  caption: [Overview of the three quantities mass, momentum and energy over time for varying grid sizes, $N$ with $alpha=50$.]
) <fig:2_d_conservation>

For small grid sizes, we see that the quantities are clearly not conserved for any value of $c$, which we attribute the phenomenon of _aliasing_, which is discussed further below in @sec:2.e. However, when we move up in grid size, we see that the estimated quantities flatten and become constants. This is what we expect from a solution to the KdV solution, and this further validates the convergence of our solver. 

== Aliasing Errors <sec:2.e>
For non-linear terms in a modelled differential equation, here taken to be the product $w(x)$ of two functions $u(x), v(x)$, we get products in the nodal domain:
$
w(x) = u(x) v(x)
$

Each of the functions $u(x), v(x)$ may be considered as an expansion in the modal basis, here taken to be the trigonometric polynomials:
$
u(x) &= sum_(k=-N/2)^(N/2-1) hat(u)_k e^(i k x), quad
v(x) &= sum_(k=-N/2)^(N/2-1) hat(v)_k e^(i k x), quad
w(x) &= sum_(k=-N/2)^(N/2-1) hat(w)_k e^(i k x)\
$

Where the coefficients for $β ∈ {u, v, w}$ are given as:
$
tilde(β)_k = 1/N sum_(j=0)^(N-1) β_j e^(-i k x_j)
$

With which we consider the modal coefficients of the product:
$
w_k
&= 1/N sum_(j=0)^(N-1) ω_j e^(-i k x_j)
&= sum_(m=-N/2)^(N/2-1) sum_(l=-N/2)^(N/2-1) hat(v)_m hat(u)_l overbrace(1/N sum_(j=0)^(N-1) e^(i (m + l -k) x_j), δ_(m+l,k+n N)\, quad n∈{-1, 0, 1})
&= hat(ω)_k + underbrace(sum_(m+l=±N) hat(v)_m hat(u)_l, "aliasing errors")
$

Which tersely outlines the treatment in @L7_slides[S.14-15] and reveals by the orthogonality condition that modes at $k±N$ are indistinguishable from the $k$th mode, resulting in _aliasing_.

As such, the element-wise multiplication in the nodal space becomes a convolutional sum in the modal space, which by construction will lead to the creation of higher-order modes, which are _aliased_ when restricted to the original lower-order discretisation. This aliasing effect pollutes the low frequency modes. 

In order to visualise this effect we investigate the the temporal evolution of the spectral coefficients, as shown in @fig:2e. By inspection of the problem, which features a travelling wave, we would expect the frequency spectrum to remain constant over time due to the constant shape of the wave – the movement of the wave is represented by the argument of the spectral coefficients, but their amplitude should remain constant.

#figure(
  image("output/2_e.png"),
  caption: [
    Temporal evolution of Fourier coefficients for the KdV equation with parameters $x ∈ [-30, 30], c=1, t ∈ [0, 3]$ with $N_"grid"=28$. Note that the green lines denote the difference between the initial and final spectra and are to be read off of the right-hand ordinate.
  ]
) <fig:2e>

Inspection of the left panel in @fig:2e shows that without any mitigation, we _create_ energy at the higher order modes as we evolve the system in time, which leads to significant aliasing errors, especially for lower spatial discretisations. This matches the results found in @sec:2.d, as seen in @fig:2_d_conservation.

In order to mitigate against this effect, we may perform the convolutional sum in a space containing higher frequency basis functions, as described by @L7_slides[S. 11—22] where aliasing effects and _Orszag's 3/2 de-aliasing rule_, which states that the product be carried out in a spectral space of order $M ≥ 3/2 N$ in order to achieve alias-free convolution operations in the spectral space.

In practice, this is by transforming the factors $v(x), u(x)$ into the Fourier space and padding them with an additional $N/2$ higher-frequency modes the associated coefficients being zero, then inverse transforming back into the physical space to carry out the spectral convolution as an element-wise multiplication, after which the product is once more transformed into the frequency space where the padded modes are removed alongside their coefficients, which now are polluted by aliasing, after which the final inverse transformation back into the physical space is undertaken. An outline of the this method implemented in Python may be found in @app:orszag.

Using the Orszag's rule implementation, we once more simulate our system to find the
evolution shown in the right panel of @fig:2e, which now shows significantly reduced aliasing. The remaining change in the spectrum are attributed as artifacts due to  the small boundary condition violation described the sections above. Note that the scale of the right axis corresponding to the difference is different in the two panels.

=== Reanalysis of conserved quantities using Orszag's de-aliasing rule

Having now implemented Orszag's 3/2 de-aliasing rule, we repeat the experiment described in @sec:2.d to produce @fig:2e_conserved_qtys.

#figure(
  image("output/2_e_conservation.png"),
  caption: [
    Overview of conserved quantities simulated with the Orszag 3/2 de-aliasing rule implemented.
  ]
) <fig:2e_conserved_qtys>

As apparent when comparing @fig:2e_conserved_qtys which is produced with the de-aliasing rule enabled and @fig:2_d_conservation in which aliasing mitigation was undertaken, we observe significantly improved stability in the system seen from the perspective of conservation of energy, momentum, and mass.

This confirms out previous suspicion that aliasing effects were to blame for the variation in these quantities as the system was evolved.

== Collision of Two Waves <sec:2.f>

We leverage the implementation for the KdV equation outlined in @sec:2.c alongside Orszag's 2/3rds dealiasing rule rule described in @sec:2.e to model the collision of two solitions
constructed using the initial condition corresponding to the superposition of two solitons:
$
u_0(x) = 1/2 [c_1 sech^2 (1/2 sqrt(c_1) (x - x_(0,1)) + c_2 sech^2 (1/2 sqrt(c_2) (x - x_(0,2))]
$

Over the domain $x ∈ [-L_x, L_x]$ with $L_x = 75$, $(x_(0,1), c_1) = (-40, 0.5)$ and $(x_(0,2), c_2) = (-15, 0.25)$ over the time interval $t∈ [0, 120]$

We visualise the solution using a ridge plot, though notably without the smoothing kernel that is often applied in such plots, in order to obtain @fig:2f.

#figure(
  image("output/2_f.png"),
  caption: [
    Unsmoothed ridge plot of two solitons colliding using the KdV equation and a Fourier collocation method for the spatial discretisation and the RK4 time-stepping algorithm. Note that the plotted domain is slightly smaller than the simulated domain to aid visibility.
  ]
) <fig:2f>

It should be noted that an exact solution may not be produced using simply the superposition of solitons due to the non-linearity of the KdV equation, which enables interaction between the colliding solitons. If the dynamics of the system had instead been described by simple advection, the evolution of the system could have been described by simple superposition of the solitons.

It was identified that the maximum of the resulting solution do not change with time and hence the stability condition from @eq:stability_region could be calculated with eigenvalue of the linearized operator $dash(L)_N$ for $t=0$.

== Scalability Analysis <sec:2.g>

To get a better understanding the scalability of our solver implementation, we benchmark a version where the differentiation is done by carrying out a Fourier transform such that the convolutional sum of the differentiation may be undertaken in $cal(O)(N)$ time as a element-wise multiplication in the modal domain, after which it is returned to the nodal domain with an inverse Fourier transformation.

Each transformation is performed using a Fast Fourier Transform algorithm, which has $N log(N)$ scaling and as such will outperform the differentiation matrix based approach, which has $N^2$ scaling for dense matrices as is the case here.

We instrument the solver to capture the time taken across the time-stepping algorithm, which solves the spatially discretised system at each time step using a Fourier collocation method, which is then divided by to number of evaluations of the right-hand side function, which in this case is the collocation method.

This approach approach yields the scaling seen in @fig:2g_internal, which notably does _not_ follow the expected scaling. This is attributed to significant and variable overhead in the time-stepping algorithm, as well as the analysis being carried out at relatively low $N$ where linear or constant time operations may still dominate the performance.

#figure(
  image("output/2_g_internal_timing.png"),
  caption: [
    Timing of solver steps collected using _internal_ timing method.
  ]
) <fig:2g_internal>

To mitigate against this and properly benchmark the solver, the Fourier collocation solver is used to solve only the spatial problem:
$
-6 u ∂_x u - ∂_(x x x) u = 0
$

This process does not involve the time-stepping algorithm at all and is able to more directly sample the performance as it would be per time step. The problem is solved $N_"repeat"=1000$ times for each grid size $N ∈ [32, 4096] ∈ ℕ$ following $N_"warmup"=100$ solutions that are discarded to prevent burn-in effects, with the result shown in @fig:2g_external.

#figure(
  image("output/2_g_external_timing.png"),
  caption: [
    Timing of solver steps collected using _external_ timing method.
  ]
) <fig:2g_external>

We note that in both cases the matrix-based solver has a step in performance at $N=64$, which we theorise to arise from exceeding a cache-level or no longer fitting in a SIMD operation. This also highlights the importance of other operations with below quadratic scaling at low $N$, which causes the matrix-based implementation to outperform the FFT-based implementation for $N≤64$.

As we go to higher grid discretisations we observe the expected scaling characteristics where the matrix-based implementation tends to $cal(O)(N^2)$ in the asymptotic limit while the FFT-based differentiation tends to $cal(O)(N log(N))$.

As hinted at by the generally higher performance of the $N = 2^n$ sizes, there are certain discretisations that perform very favourably in FFT algorithms. A deeper discussion of these is beyond the scope of this report, but as shown in @fig:2g_external it may often be beneficial to use an $N=2^n, n ∈ ℕ_1$ discretisation to achieve the best performance.

Our implementation features standard techniques for improved performance including vectorisation of operations where possible. Notably, for the matrix-based implementation the differentiation matrices are precomputed to prevent superfluous computation at every evaluation of the spatial problem. Just-in-time optimisation using `numba` was also undertaken, but as evidenced by @fig:2g_external, we see little performance gains from this, which may be understood by realising that the vectorised `numpy` operations are highly optimised already and leave little on the table performance-wise. A last optimisation that is undertaken is preallocation of the memory for the arrays, which saves slightly on memory allocation and deallocation. This trick is usually very significant for GPU-based setups, where the memory is passed back on forth between devices, but does also provide some speedup in this case.

== Solving an IVP with Space-Time Discretisation <sec:2.h>
We consider the linear advection equation which is an Initial Value Problem
$
  pdv(phi, t) + a pdv(phi, x) &= 0, quad x in thin ]0, 2pi[,thin thin t>0 \
  phi(x, 0) &= phi_0(x), quad x in [0, 2pi] \
  phi(0, t) &= g_l (x), quad t >= 0 \  
$
We choose to consider the soliton function from @sec:2.c as the solution to the linear advection equation, restricted to the domain $[0, 2pi]$ and we further derive from it the boundary conditions $phi_0, g_l$. The function is the following
$
  phi(x, t) = f(x - c t), quad f(x)=1/2 c sech^2(1/2sqrt(c)(x-x_0)), quad x in [0, 2pi].
$
We wish to solve the Initial Value Problem as a Boundary Value Problem with the spectral Legendre Collocation Method. We consider the formulation of the linear advection equation as the Boundary Value Problem
$
  pdv(phi, t) + a pdv(phi, x) &= 0, quad (x, t) in thin ]0, 2pi[ thin times thin ]0, T_max]\
  phi(x, 0)& = phi_0(x), quad x in [0, 2pi] \
  phi(0, t)& = g_l (x), quad t in [0, T_max]
$
where we denoted by $T_max$ the temporal end of the domain.

We consider the following coordinate transformations to map the two dimensional domain of the BVP problem to the $(z_x, z_t) in [-1,1]^2$ domain
$
  t = (T_max)/2 z_t + T_max / 2, quad & pdv(, t) = 2 / (T_max) pdv(, z_t) \
  x = pi z_x + pi, quad  &pdv(, z) = 1 / (pi) pdv(, z_x).
$
After coordinate transformation we obtain the following BVP
$
  2 / (T_max) pdv(phi, z_t) (z_x, z_t)  + a 1/ (pi) pdv(phi, z_x) (z_x, z_t)& = 0, quad (z_x, z_t) in [-1,1]^2\
  phi(z_x, 0) &= phi_0(z_x), quad z_x in [-1, 1] \
  phi(0, z_t) &= g_l (z_t), quad z_t in [-1, 1].
$
We consider the modal and nodal expansion
$
  phi_N (z_x, z_t) = sum_(i=0)^N_x sum_(j=0)^N_t hat(phi)_(i j) P_i (z_x) P_j (z_t) = sum_(i=0)^N_x sum_(j=0)^N_t phi_(i j) h_i (z_x) h_j (z_t)
$
where the grid points are the Jacobi Gauss-Lobatto grid points, and the $P$ are the Legendre polynomials. We use the spectral Legendre Collocation Method, and we therefore want the residual to be zero at each collocation node,
$
  lr({ 2 / (T_max) (phi_N)_(z_t) + a / (pi) (phi_N)_(z_x)}  |)_((z_x^i, z_t^j)) = 0, wider i=0,1,...,N_x, quad j = 0,1,...,N_t
$
we furhter impose appropriate boundary conditions. We construct the partial derivatives of $phi_N$ by using the Differentiation matrices
$
  pdv(phi_N, z_x)(z_x, z_t) = D_(z_x) bold(phi), quad pdv(phi_N, z_t)(z_x, z_t) = D_(z_t) bold(phi), quad D_(z_x) in bb(R)^(N_x times N_x), quad D_(z_t) in bb(R)^(N_t times N_t)
$
We construct the global operator by using the Kronecker products with the Identity matrices
$
  L_N = ((2 / (T_max) D_(z_t)) times.o I_(z_x)) + (I_(z_t) times.o (a / (pi) D_(z_x)))
$
after incorporating the boundary conditions the system matrix becomes $tilde(L)_N$, and the right hand side of the system equation becomes $tilde(f)$. The system of equations which yields the nodal solution $phi_N$ at collocation points becomes
$
  tilde(L)_N phi_N = tilde(f)
$
We solve the problem with the following parameters $T_max = 1, c = 10, a = 1, x_0 = pi$. The convergence plots are presented in @fig:2h_convergence. 

We note the relationship between the temporal discretisation $N_t$ and the spatial discretisation $N_x$. Based on the plots it could be estimated that we should ensure $N_x > gamma N_t$ where $gamma approx 4$, when the gird is refined. It was verified that $gamma$ does not depend on the problem parameters. Consider the convergence curve for fixed $N_t$ and as a function of $N_x$. When $N_x < gamma N_t$ the spatial error dominates the overall error, since the temporal dimension have been discretised enough, and we get the spectral convergence based on the semilog plot. Once $N_x > gamma N_t$, further spatial refinement do not result in smaller norm of the error, which indicates that the temporal error dominates. This statement is further supported by the fact that with higher fixed $N_t$ the plateau in the convergence plot is seen on lower values of the norm of the error. 

The snapshots of the solution are presented in @fig:2h_sample_solution. 

#figure(
  image("output/2_h_convergence.png"),
  caption: [
    Convergence plot for the linear advection equation solved as the BVP with the Spectral Legendre Collocation method.
  ]
) <fig:2h_convergence>

#figure(
  image("output/2_h_advection_solution.png"),
  caption: [
    Snapshots of the solutions for the linear advection equation in the original domain $[0, 2pi]$.
  ]
) <fig:2h_sample_solution>

// General coordinate transform
// $
//   x = (b - a)/2 z + (b + a)/2
// $
// Consider
// $
  
// $
// apply boundaries for $[-infinity, infinity]$, to be $[-alpha, alpha]$.
// $
//   pdv(, x) = 1 / (alpha) pdv(, z_x)
// $
// we get
// $
//   2 / (T_max) pdv(phi, z_t) + a 1/ (alpha) pdv(phi, z_x) = 0
// $

= Appendix

== Outline of Orszag's 3/2 de-aliasing rule <app:orszag>

```python
# Compute wave numbers and differentiation operators
k = 2 * np.pi * np.fft.fftfreq(N_grid, d=(L / N_grid))  # Wave numbers
D_hat = 1j * k
D3_hat = (1j * k) ** 3

# Pad
N_size = u_hat.size
M_size = int(np.ceil(N_size * 3 / 2))  # Pad size
if N_size % 2 == 0:
    N_pos = N_size // 2 + 1  # Includes 0, positive freqs, Nyquist freq
    N_neg = (
        N_size // 2 - 1
    )  # Includes negative freqs excluding Nyquist freq
else:
    N_pos = (N_size + 1) // 2  # Includes 0, positive freqs
    N_neg = N_pos - 1  # Includes negative freqs

u_hat_padded = np.zeros(M_size, dtype=np.complex128)

# Add in original coefficients, note new higher frequencies are at center
u_hat_padded[:N_pos] = u_hat[:N_pos]
u_hat_padded[M_size - N_neg :] = u_hat[N_size - N_neg :]

# Differentiate in spectral space
k_padded = 2 * np.pi * np.fft.fftfreq(M_size, d=(L / M_size))
D_hat_padded = 1j * k_padded
u_x_hat_padded = D_hat_padded * u_hat_padded

# Transform back to physical space
u_padded = np.fft.ifft(u_hat_padded)
u_x_padded = np.fft.ifft(u_x_hat_padded)

# Compute product in physical space
uu_x_padded = u_padded * u_x_padded

# Now remove padding in spectral space
uu_x_hat_padded = np.fft.fft(uu_x_padded)
uu_x_hat = np.zeros(N_size, dtype=np.complex128)
uu_x_hat[:N_pos] = uu_x_hat_padded[:N_pos]
uu_x_hat[N_size - N_neg :] = uu_x_hat_padded[M_size - N_neg :]

# Fix normalisation gain made by FFT computation
uu_x_hat *= M_size / N_size  # Should be ≈ 3/2

u_xxx_hat = D3_hat * u_hat

u_alias_free = np.fft.ifft(-6 * uu_x_hat - u_xxx_hat)
```

#bibliography("references.bib")
