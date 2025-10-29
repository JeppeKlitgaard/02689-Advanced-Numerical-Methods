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

In either case, we consider the coordinate transformation at hand:
$
x = (b - a)/2 z + (a+b)/2
$

With associated differentials:
$
dif x = (b - a)/2 dif z
$

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

Which gives modal and nodal solutions:
$
hat(vv(u)) &= (cal(L)_N)^(-1) hat(vv(f))\
u (x) &= sum_(n=0)^N hat(u)_n ϕ_n (x)
$

In order to find the discretisation given by @eq:1_lin_sys_1, we may obtain a sparse matrix $mm(A)$ by finding a recurrence relation dictating the coefficients of the algebraic equations

TODO MERGE WITH BELOW

---

note that we need to apply coordinate transformation but it cancels out on both side of the weak form

weight function for Legendre w(x) = 1, so convenitent

Consider
$
  a hat(u)^((2))_n + b hat(u)^((1))_n + c hat(u)_n= hat(f)_n; quad hat(u)^((2))_n =
  -(b)/(a) hat(u)^((1))_n - c/a hat(u)_n + 1/a hat(f)_n
$
Recursion
$
  hat(u)^((q-1))_n = 1/(2n-1) hat(u)^((q))_(n-1) - 1/(2n+3) hat(u)^((q))_(n+1)
$
Consider $q=2, q=1$
$
  hat(u)^((1))_n = 1/(2n-1) hat(u)^((2))_(n-1) - 1/(2n+3) hat(u)^((2))_(n+1); quad hat(u)_n = 1/(2n-1) hat(u)^((1))_(n-1) - 1/(2n+3) hat(u)^((1))_(n+1)
$
Plug in
$
  hat(u)^((1))_n = 1/(2n-1) (-(b)/(a) hat(u)^((1))_(n-1) - c/a hat(u)_(n-1) + 1/a hat(f)_(n-1)) - 1/(2n+3) (-(b)/(a) hat(u)^((1))_(n+1) - c/a hat(u)_(n+1) + 1/a hat(f)_(n+1))
$
Rewrite
$
    hat(u)^((1))_n = -(b)/(a) (1/(2n-1) hat(u)^((1))_(n-1) - 1/(2n+3) hat(u)^((1))_(n+1) ) - c/a (1/(2n-1) hat(u)_(n-1) - 1/(2n+3) hat(u)_(n+1) ) + 1/a (1/(2n-1) hat(f)_(n-1) - 1/(2n+3) hat(f)_(n+1))
$
Use recursion
$
    hat(u)^((1))_n = -(b)/(a) hat(u)_n - c/a (1/(2n-1) hat(u)_(n-1) - 1/(2n+3) hat(u)_(n+1) ) + 1/a (1/(2n-1) hat(f)_(n-1) - 1/(2n+3) hat(f)_(n+1))
$
Consider
$
  1/(2n-1) hat(u)^((1))_(n-1) &= 1/(2n-1) (-(b)/(a) hat(u)_(n-1) - c/a (1/(2n-3) hat(u)_(n-2) - 1/(2n+1) hat(u)_(n) ) + 1/a (1/(2n-3) hat(f)_(n-2) - 1/(2n+1) hat(f)_(n))) = \
  &= (- c/a 1/(2n-1) 1/(2n-3)) hat(u)_(n-2) + (- b/a 1/(2n-1)) hat(u)_(n-1) + (c/a 1/(2n-1) 1/(2n+1)) hat(u)_(n) \
  &+ (1/a 1/(2n-1) 1/(2n-3)) hat(f)_(n-2) + (- 1/a 1/(2n-1) 1/(2n+1)) hat(f)_(n)

$
Consider
$
  1/(2n+3) hat(u)^((1))_(n+1) &= 1/(2n+3) (-(b)/(a) hat(u)_(n+1) - c/a (1/(2n+1) hat(u)_(n) - 1/(2n+5) hat(u)_(n+2) ) + 1/a (1/(2n+1) hat(f)_(n) - 1/(2n+5) hat(f)_(n+2))) = \
  &= (- c/a 1/(2n+3) 1/(2n+1)) hat(u)_(n) + (- b/a 1/(2n+3)) hat(u)_(n+1) + (c/a 1/(2n+3) 1/(2n+5)) hat(u)_(n+2) \
  &+ (1/a 1/(2n+3) 1/(2n+1)) hat(f)_(n) + (- 1/a 1/(2n+3) 1/(2n+5)) hat(f)_(n+2)
$
The terms which are overlaping are $hat(u)_(n), hat(f)_(n)$. Consider
$
  (c/a 1/(2n-1) 1/(2n+1)) hat(u)_(n) - (- c/a 1/(2n+3) 1/(2n+1)) hat(u)_(n) - hat(u)_(n) = (c/a 1/(2n-1) 1/(2n+1) + c/a 1/(2n+3)1/(2n+1) - 1) hat(u)_(n) =\
  ( c/a 2 / ((2n-1)(2n+3)) - 1)hat(u)_(n)
$
Consider
$
  (- 1/a 1/(2n-1) 1/(2n+1)) hat(f)_(n) - (1/a 1/(2n+3) 1/(2n+1)) hat(f)_(n) = - (1/a 1/(2n-1) 1/(2n+1) + 1/a 1/(2n+3) 1/(2n+1)) hat(f)_(n) = \
  - 1/a 2 / ((2n-1)(2n+3)) hat(f)_(n)
$

Finally
$
  1/(2n-1) hat(u)^((1))_(n-1) - 1/(2n+3) hat(u)^((1))_(n+1) - hat(u)_(n) = 0
$
We get
$
  a_(n, n-2) hat(u)_(n-2) + a_(n, n-1) hat(u)_(n-1) + a_(n, n) hat(u)_(n) + a_(n, n+1) hat(u)_(n+1) + a_(n, n+2) hat(u)_(n+2) = g_(n, n-2) hat(f)_(n-2) + g_(n, n) hat(f)_(n) + g_(n, n+2) hat(f)_(n+2)
$
The coefficients are the following
$
  a_(n, n-2) &= - c/a 1/(2n-1) 1/(2n-3) \
  a_(n, n-1) &= - b/a 1/(2n-1) \
  a_(n, n)   &= c/a 2 / ((2n-1)(2n+3)) - 1 \
  a_(n, n+1) &= b/a 1/(2n+3)\
  a_(n, n+2) &= - c/a 1/(2n+3) 1/(2n+5)\
  g_(n, n-2) &= - 1/a 1/(2n-1) 1/(2n-3)\
  g_(n, n) &= 1/a 2 / ((2n-1)(2n+3))\
  g_(n, n+2) &= - 1/a 1/(2n+3) 1/(2n+5)\

$

OMG are we supposed to do the coordinate transformation here ... :(

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
Where we have used $w(x) = (1-x)^α (1+x)^β = 1$ for the Legendre polynomials characterised by $α=β=0$.

This yields the collocation condition:
$
u(x_i) = sum_(n=0)^N hat(u)_n ϕ_n (x_i) wider i = 0, 1, ..., N
$ <eq:1_collocation>

Where we may determine the expansion coefficients ${hat(u)_n}_(n=0)^N$ by orthogonality:
$
hat(u)_N = ((u, ϕ_n)_(w, N))/((ϕ_n, ϕ_n)_(w, N)) = 1/γ_n ∫_a^b u_N (x) ϕ_n (x) dif x wider n = 0, 1, ..., N
$ <eq:1_modal_expansion_coeffs>

Where the integral will usually be computed using a quadrature. If we are free to pick our collocation points, $x_i$, and quadrature rule, we may pick the Legendre-Gauss-Lobatto points and the Gaussian quadrature a result that is exact for orthogonal polynomials up to degree $2n-1$.

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

=== Discussion of Methods <sec:1.a_discussion>


== Irrotational flow around a cylinder <sec:1.b>

We are given the Laplace's Equation in 2 dimensions expressed in polar coordinates:
$
∇^2 ϕ = 1/r pdv(, r) (r pdv(ϕ, r)) + 1/r^2 pdv(ϕ, θ, 2) = 0
$

---

Consider
$
  1/r partial / (partial r) (r (partial phi)/ (partial r) (r, theta)) + 1/r^2 (partial^2 phi) / (partial theta ^ 2) (r, theta) = 0 \
  r (partial phi) / (partial r) (r, theta) + r^2 (partial^2 phi)/ (partial r^2) (r, theta) + (partial^2 phi) / (partial theta ^ 2) (r, theta) = 0 \
  phi(r_1, theta) = g_1(theta); quad phi(r_2, theta) = g_2(theta) \
  phi(r, 0) = phi(r, 2 pi) = g_3(r)
$
Exact solution
$
  phi(r, theta) = V_infinity (r + r_1^2 / r) cos(theta)
$
Consider
$
  r = (r_2 - r_1)/2 z_r + (r_2 + r_1)/2, quad partial / (partial r) = 2/(r_2 - r_1) partial / (partial z_r), quad partial^2 / (partial r^2) = 4/(r_2 - r_1)^2 partial / (partial z_r)
$
Define
$
  alpha_1 = (r_2 - r_1)/2, alpha_2 = (r_2 + r_1)/2, quad r = alpha_1 z_r + alpha_2, quad partial / (partial r) = 1/alpha_1 partial / (partial z_r), quad partial^2 / (partial r^2) = 1/alpha_1^2 partial / (partial z_r)
$
Consider
$
  theta = (2 pi)/2 z_theta + (2 pi)/2 = pi z_theta + pi, quad partial / (partial theta) = 1/pi partial / (partial z_theta), quad partial^2 / (partial theta^2) = 1/pi^2 partial^2 / (partial z_theta^2)
$
Consider
$
    r (partial phi) / (partial r) + r^2 (partial^2 phi)/ (partial r^2) + (partial^2 phi) / (partial theta ^ 2) = 0 \

    (alpha_1 z_r + alpha_2) 1/alpha_1 partial / (partial z_r) phi + (alpha_1 z_r + alpha_2)^2 1/alpha_1^2 partial^2 / (partial z_r^2) phi + 1/pi^2 partial^2 / (partial z_theta^2) phi = 0 \

      (z_r + alpha_2 / alpha_1) partial / (partial z_r) phi + (z_r + alpha_2 / alpha_1)^2 partial^2 / (partial z_r^2) phi + 1/pi^2 partial^2 / (partial z_theta^2) phi = 0
$
Solution becomes
$
  phi(z_r, z_theta) = (alpha_1 z_r + alpha_2 + r_1^2 / (alpha_1 z_r + alpha_2)) cos(pi z_theta + pi)
$
Collocation method discretization
$
  D_x = "diag"(beta_1, ..., beta_N) D_r + "diag"(beta_1^2, ..., beta_N^2) D_r^2\
  D_y = 1/pi^2 D_theta^2\
  L_N = "kron"(I, D_x) + "kron"(D_y, I)
$
We define the new coordinate system $(z_r, z_theta) = (x, y)$, we use the row-wise ordering.

Two dimensional Vandermore matrices
$
  U = V_x hat(U) V_y^T
$

= Time-Dependent Problems <sec:2>

#counter(heading).update((2, 2))
== <sec:2.c>
== <sec:2.d>
== <sec:2.e>
== <sec:2.f>
== <sec:2.g>
== <sec:2.h>
