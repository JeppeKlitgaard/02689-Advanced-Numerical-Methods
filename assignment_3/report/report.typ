#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "@preview/cheq:0.2.2": checklist
#import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1"
#import "@preview/algorithmic:0.1.0"
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT, pdv
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

#set document(title: "Assignment 3")
#align(center, [
  #title()

  #set text(size: 14pt)
  02689 — Advanced Numerical Methods for Differential Equations
])

#v(8pt)
- Jeppe Klitgaard <`s250250@dtu.dk`>
- Tymoteusz Barcinski <`s221937@dtu.dk`>
- Pernille Christie <`s204249@dtu.dk`>
/*
= TODO

- Energy calculations
- Upwind scheme
- Skew-symmetric form
- Maybe also do square domain?
- Poster

= Problem Formulation <sec:problem>

- Derive a SEM solver for the rotating hill problem. SEM requires a special type of discretization that is skew-symmetric for stability. This form is a result of a rewrite of the mathematical model using vector calculus, and then it is discretized.
- Solve the problem with a SEM in a circular domain (curved boundaries), i.e. a high-order scheme that has high-order and low-order in the same formulation.
- Perform convergence tests, i.e. h-convergence and p-convergence tests. 
- Use the solver to demonstrate break-even points for when high-order is more cost-efficient than low-order vs. high-order as a function of temporal integration time using an explicit 4th order runge-kutta method.
- Derive skew-symmetric form
- Prove energy conservation

Questions for Allan:
- Governing equation (is there a document or paper we should use as reference?)
- Could he elaborate a bit on first point, hints on rewrite of mathematical model to get skew symmetry
- Hints/tools to do the mesh for a circular domain
  - Arises from the fact that rotation matrices are skew symmetric?
- What exactly is the difference between continuous and discontinuous Galerkin?
  - Intuition: how can we get $K N^p$ scaling if we impose continuous BC's between elements? Does this not couple the points of elements together?
- Should we use continuous or discontinuous method?

Hints:
- Upwind boundary conditions
- Create grid for each element
- Global indices for nodes in grid

For week 12:
- Connectivity grid
- Global index

Triangles - page 171, chapter 6.1*/

= Skew-Symmetric Split Form of the Advection Equation
We consider the conservative advection equation
$
pdv(u, t) + nabla dot (vv(a)u)=&0, quad &(x,y)in Omega,thin t>0\
u(x,y,0) =& u_0(x,y),  quad & (x,y)in Omega
$ <eq:standard_eq>
where $vv(a)$ is some velocity field. We assume that the flow-field is divergence free, i.e.
$nabla dot vv(a) = 0.$

/*
where
$
vv(a) = (2pi y_c, -2pi x_c),
$
with 
$
x_c (t)&=x_0 cos(2 pi t) + y_0 sin(2pi t),\
y_c (t)&=y_0 cos(2pi t) - x_0 sin(2pi t).
$*/

We use the chain rule to rewrite @eq:standard_eq,
$
pdv(u,t)+ nabla dot (vv(a)u) = pdv(u,t)+vv(a) dot nabla u+(nabla dot vv(a))u=pdv(u,t)+vv(a) dot nabla u = 0.
$
By taking the average of the two equations, we get
$
pdv(u,t)+1/2[nabla dot (vv(a) u)+vv(a) dot nabla u]=0.
$ <eq:skew_symmetric>
We show below that this results in a skew-symmetric operator.

= Weak Formulations: Conservative and Split Form

== Separate Volume and Boundary Contribution
To derive the weak form of @eq:standard_eq, we multiply by some test function $v$ and take the integral over the domain $Omega$,
$
integral_(Omega)v pdv(u,t) dif Omega +integral_Omega v nabla dot (vv(a) u) dif Omega=0.
$
We can then use integration by parts to rewrite the second term, yielding
// $
// integral_(Omega)v pdv(u,t) dif Omega +integral_(partial Omega) v u(vv(a) dot vv(n) )dif s - integral_Omega (vv(a)u) dot nabla v dif Omega=0.
// $
// JK: I reordered some terms and factors to get it to match (11) in assignment description
$
integral_(Omega)v pdv(u,t) dif Omega - integral_Omega nabla v dot (vv(a)u) dif Omega + integral_(partial Omega) v (vv(a) dot vv(n)) u dif s =0.
$ <eq:weak_form_normal>

== Separate Volume and Boundary Contribution for Split Form
Doing the same thing for @eq:skew_symmetric, we get 
$
integral_(Omega)v pdv(u,t) dif Omega +1/2[integral_Omega v nabla dot (vv(a) u) dif Omega + integral_Omega v vv(a) dot nabla u dif Omega]=0.
$
Again using integration by parts on the first term of the parenthesis, we get
$
integral_(Omega)v pdv(u,t) dif Omega + 1/2[integral_(partial Omega) v u(vv(a) dot vv(n) )dif s -integral_Omega u vv(a) dot nabla v dif Omega + integral_Omega v vv(a) dot nabla u dif Omega]=0.
$ <eq:weak_form_split>

// PC: Still very much in doubt of whether to use this or not - right now I am leaning to not
/*If we choose basis functions that vanish at the boundary, the boundary integrals will be 0, so the final expression will be
$
integral_(Omega)v pdv(u,t) dif Omega + 1/2[integral_Omega v vv(a) dot nabla u dif Omega - integral_Omega u vv(a) dot nabla v dif Omega]=0.
$ <eq:skew_weak>*/


// JK: I think that the basis functions don't 'disappear' in the sense that they evaluate to zero, but more in the sense that they evaluate to exactly the boundary condition at the boundary. Since we have Lagrange basis, they will become exactly the imposed boundary condition at the boundaries and will only disappear if the boundary conditions mean that the integral ends up as zero.
// PC: I have based this on the derivation on slide 13 in lecture 10

== Decomposition of Domain
To ease computations, we define a mesh which divides $Omega$ into multiple smaller domains, $Omega_k$, such that 
$
Omega approx tilde(Omega) = union.big_(k=1)^K Omega_k.
$
Since the integral operator is linear, we can decompose each area integral into a sum of integrals over the $Omega_k$'s. This yields
$
sum_(k=1)^K [ integral_(Omega_k)v pdv(u,t) dif Omega_k + 1/2(integral_Omega_k v vv(a) dot nabla u dif Omega_k - integral_Omega_k u vv(a) dot nabla v dif Omega_k)] + 1/2integral_(partial tilde(Omega)) v u(vv(a) dot vv(n) )dif s=0.
$ 
// PC: Unsure of how to write this if we have the boundary integrals
We will denote the domains with edges along the boundary by $tilde(Gamma)$. We now require the residual to vanish in each interior domain, that is,
$
 integral_(Omega_k)v pdv(u,t) dif Omega_k + 1/2(integral_Omega_k v vv(a) dot nabla u dif Omega_k - integral_Omega_k u vv(a) dot nabla v dif Omega_k)=0, quad forall Omega_k in.not tilde(Gamma).
$  <eq:system_omegak>
For the boundary elements of the mesh, we will include the boundary terms  
$
integral_(Omega_j)v pdv(u,t) dif Omega_j + 1/2(integral_Omega_j v vv(a) dot nabla u dif Omega_j - integral_Omega_j u vv(a) dot nabla v dif Omega_j + integral_(partial Omega_j inter partial tilde(Omega)) v u(vv(a) dot vv(n) )dif s), quad forall Omega_j in tilde(Gamma)
$
which we will consider in the following sections.

= Upwind Numerical Flux at the Boundary

We describe the numerical flux at the boundary $∂Ω$ by
$
f_n = (vv(a) ⋅ vv(n)) u .
$

We let $u^-$ and $u^"bc"$ denote the value of $u$ as given by the interior trace or the prescribed boundary condition respectively, which leads to the following definition of the _upwind scheme_:
$
  f_n^* (u^-, u^"bc") = 1/2 (a_n + |a_n|) u^- + 1/2 (a_n - |a_n|) u^"bc" .
$ <eq:3_flux>

== Explanation of Notation

As such, the flux $f^*_n$ will be determined by $u^-$ or $u^"bc"$ depending on the advection speed $a_n$ along the normal $vv(n)$, as shown in @fig:3_a.

#figure(
  image("3_a.jpg", height: 4cm),
  caption: [
    Illustration of geometrical understanding of $u^-$, $u^"bc"$, and $a_n$. Outward flow is shown.
  ]
) <fig:3_a>

== Inflow and Outflow Fluxes

As seen in @fig:3_a, the flux across the boundary $∂Ω$ can either be _inward_ or _outward_, depending on the sign of $a_n$, where an _inflow_ flux corresponds to $a_n<0 ⇔ abs(a_n) = -a_n$ and an _outflow_ flux corresponds to $a_n>0 ⇔ abs(a_n) = a_n$ as shown.

We evaluate the fluxes for both cases:
$
f_(n,"outflow")^*(u^-, u^"bc") &= 1/2 (a_n + a_n) u^- + cancel(1/2 (a_n - a_n) u^"bc") = a_n u^-,
wider & a_n > 0
\
f_(n,"inflow")^*(u^-, u^"bc") &= cancel(1/2 (a_n - a_n) u^-) + 1/2 (a_n - (-a_n)) u^"bc" = a_n u^"bc",
wider & a_n < 0
$ <eq:upwind_flux>

== Weak Form and Boundary Condition Imposition

We may now rewrite the boundary terms of @eq:weak_form_normal and @eq:weak_form_split using this upwind flux:
$
∫_(∂Ω) v (vv(a) ⋅ vv(n)) u dif s = ∫_(∂Ω) v f_n dif s = ∫_(∂Ω) v f_n^* dif s 
$

This is helpful in the imposition of boundary conditions as it allows us to determine which boundary nodes should be determined using the boundary condition $u^"bc"$, which by inspection of @eq:upwind_flux becomes those satisfying the _inflow condition_, $a_n < 0$.

// NOTE, JK: I don't understand why this is helpful. Rotational symmetry implies $vv(a) ⋅ vv(n) = 0$ everywhere on $∂Ω$? We should return to this section and write a bit more once we figure out if this is all just busywork to conclude that the boundary flux must be zero, or whether this bit is actually significant.

= Semi-Discrete Spectral Element Formulation and Global Matrices

== Transforming the Domain
As mentioned above, we divide $Omega$ into a mesh of multiple triangles. These triangles must be transformed into a reference domain, $I={vv(r)=(r,s) | (r,s)>=-1, s+r <= 0}$, which we can perform calculations on, see // PC: insert figure of transformation
/*
#figure(
  polygon(
  fill: blue.lighten(80%),
  stroke: blue,
  (10%, 0pt),
  (30%, 2cm),
  (0%,  2cm),
)
)*/

The mapping between a point on the general triangle $vv(x)in D$ to the corresponding point on the reference triangle $vv(r) in I$ is
$
vv(Psi)(r,s)=1/2(-(r+s)vv(v)^1+(r+1)vv(v)^2+(s+1)vv(v)^3).
$
We are able to define the discrete operations on element $I$ as follows
$
cal(D)_r= cal(V)^(-1)cal(V)_r, quad cal(D)_s=cal(V)^(-1)cal(V)_s,\
vv(x)_r=cal(D)_r vv(x), quad vv(x)_s=cal(D)_s vv(x),\
vv(y)_r=cal(D)_r vv(y), quad vv(y)_s=cal(D)_s vv(y),\
vv(J)=vv(x)_r times.square vv(y)_s - vv(x)_s times.square vv(y)_r.
$
We can use these to determine the discrete directional mappings needed for later computations,
$
cal(D)_x = "diag"(vv(r)_x)cal(D)_r+"diag"(vv(s)_x)cal(D)_s,\
cal(D)_y = "diag"(vv(r)_y)cal(D)_r+"diag"(vv(s)_y)cal(D)_s,
$
where $vv(r)_x = vv(y)_s \/ vv(J)$, $vv(r)_y = -vv(x)_s \/ vv(J)$, $vv(s)_x =- vv(y)_r \/ vv(J)$ and $vv(s)_y = vv(x)_r \/ vv(J)$. 

Now for the mass-matrix, we consider the classic integral over $D$ and transform it
$
integral_D phi.alt_n (vv(x)) phi.alt_m (vv(x)) dif vv(x)=integral_I phi.alt_n (vv(r)) phi.alt_m (vv(r)) J dif vv(r).
$
Thus, for each element $Omega_k$, we will define 
$
cal(M)^((k)) = "diag"(vv(J)^((k))) cal(M),
$
where $cal(M)$ is the mass matrix for $I$ and $vv(J)^((k))$ is the Jacobian for the transformation from $Omega_k$ to $I$. 

== Semi-discrete System on Element-level
We will divide $u$ into several elements $u^((k))$ such that $u(t, vv(x))=plus.o.big_k u^((k))(t, vv(x))$, with
$
u^((k))(t,vv(x))=cases(
  u(t, vv(x))\, quad &vv(x)in Omega_k,
  0\, &"otherwise"
).
$
This means that we can restate @eq:system_omegak as
$
 integral_(Omega_k)v pdv(u^((k)),t) dif Omega_k + 1/2(integral_(Omega_k) v vv(a) dot nabla u^((k)) dif Omega_k -integral_(Omega_k) u^((k)) vv(a) dot nabla v dif Omega_k)=0, quad k = 1,...,k.
$
For each element, we will approximate the function as 
$
u^((k))(t,vv(x)) approx u_h^((k))(t, vv(x))=sum_(n=1)^(N_p) u^((k))_n (t) phi.alt_n (vv(x)).
$
For this project, we will use a Galerkin method, meaning that our test functions are the same as our basis functions $v=phi.alt_m$. Using this, we will rewrite the terms of the weak formulation. Starting with the time-derivative term,
$
(pdv(u^((k))_h, t), phi.alt_m)_(N_p)=&integral_(Omega_k)phi.alt_m (vv(x)) pdv( ,t) (sum_(n=1)^(N_p) u^((k))_n (t) phi.alt_n (vv(x))) dif vv(x)\
=&sum_(n=1)^(N_p) dv(u^((k))_n ,t) integral_(Omega_k)phi.alt_n (vv(x))phi.alt_m (vv(x)) dif vv(x)\
//=&sum_(n=1)^(N_p) dv(u^((k))_n ,t) cal(J)^((k)) integral_(I)phi.alt_n (vv(r))phi.alt_m (vv(r)) dif vv(r)\
=&sum_(n=1)^(N_p) dv(u^((k))_n ,t) cal(M)^((k))_(m n) \
=&[cal(M)^((k)) dv(vv(u)^((k)),t)]_m.
$
For the other non-boundary terms, we get
$
(vv(a) dot nabla u_h^((k)), phi.alt_m)_(N_p)=&integral_(Omega_k) phi.alt_m (vv(x)) vv(a) dot nabla (sum_(n=1)^(N_p) u^((k))_n (t) phi.alt_n (vv(x)))   dif vv(x)\
=& sum_(n=1)^(N_p) u^((k))_n (t) integral_(Omega_k)  phi.alt_n (vv(x)) (a_x pdv(,x) phi.alt_n (vv(x))+ a_y (t)pdv(,y) phi.alt_n (vv(x))) dif vv(x)\
=& [("diag"(a_x) cal(M)^((k)) cal(D)_x + "diag"(a_y) cal(M)^((k)) cal(D)_y) vv(u)^(k)]_m,
$
and
$
(u_h^((k)), vv(a) dot nabla phi.alt_m)_(N_p)=&integral_(Omega_k) (sum_(n=1)^(N_p) u^((k))_n (t) phi.alt_n (vv(x))) (vv(a) dot nabla phi.alt_m (vv(x)))  dif vv(x)\
=& sum_(n=1)^(N_p) u^((k))_n (t) integral_(Omega_k) phi.alt_n (vv(x)) (a_x lr(pdv(phi.alt_m, x)|)_(vv(x))+ a_y (t) lr(pdv(phi.alt_m, y)|)_(vv(x))) dif vv(x)\
//=&2pi sum_(n=1)^(N_p) u^((k))_n (t) ( y_c (t) cal(J)^((k)) integral_(I) phi.alt_n (vv(r)) lr(pdv(phi.alt_m, x)|)_(vv(r)) dif vv(r) - x_c (t)  cal(J)^((k)) integral_(I) phi.alt_n (vv(r)) lr(pdv(phi.alt_m, y)|)_(vv(r)) dif vv(r))\
=& [("diag"(a_x) cal(D)_x^"T"cal(M)^((k)) + "diag"(a_y) cal(D)_y^"T"cal(M)^((k))) vv(u)^((k))]_m\
=& [("diag"(a_x) (cal(M)^((k)) cal(D)_x)^"T" + "diag"(a_y) (cal(M)^((k))cal(D)_y)^"T") vv(u)^((k))]_m.
$


//AT LEAST ACCORDING TO SLIDES THE ORDER OF D_X AND M SHOULD BE SWITCHED?
// PC: You are correct, mysterious writer, and I have fixed it

Combining these expressions, we get the following system for each interior element, $Omega_k$,
$
cal(M)^((k)) dv(vv(u)^((k)),t)=&cal(L)^((k))vv(u)^((k)),
$
where
$
cal(L)^((k)):=-1/2 [thin "diag"(a_x) (cal(M)^((k))cal(D)_x - (cal(M)^((k))cal(D)_x)^"T")  + "diag"(a_y) (cal(M)^((k))cal(D)_y - (cal(M)^((k))cal(D)_y)^"T") ]
$
Note that $cal(L)^((k))$ is skew-symmetric, that is $cal(L)^((k)) = -(cal(L)^((k)))^"T"$.

Furthermore, we need to impose the weak boundary conditions, i.e. the boundary integrals. We will calculate these using Gauss-Lobatto-Legendre (GLL) quadrature,
$
1/2integral_(partial Omega_j inter partial tilde(Omega)) v u(vv(a) dot vv(n) )dif s =1/2 integral_(partial Omega_j inter partial tilde(Omega)) v f_n^* dif s =1/2 integral_(-1)^1 v f_n^* J^((j)) dif xi approx 1/2 sum_(i=1)^N_p w_i f^*_n J^((j))=:tilde(f)^((j)),
$
where we use $J^((j))$ to denote the Jacobian for the transformation from the variable running along the edge to $xi in [-1,1]$ which is the interval where Legendre polynomials are defined. For each node on $partial Omega_j inter partial tilde(Omega)$, we add this to the right hand side of the equation. We will explain exactly how this is done globally in the next section.

== Global assembly
For the global assembly, we assign each node in each element a unique global index. We can then assemble the global system and mass matrices as follows
$
cal(M)_(i j) = cal(M)_(ell n)^((k))
$
where $(i,j)$ is the global index of the vertex with local index $(ell,n)$ in element $k$. We do the same thing for the global $cal(L)$ matrix. 

To impose the week boundary conditions we need to add the $tilde(f)^((j))$'s appropriately. To do this we create a zero-vector of length $N_p K$, $tilde(vv(f))$, and then for each boundary element $Omega_j in tilde(Gamma)$, we add $tilde(f)^((j))$ to the global index of each of the boundary nodes from element $Omega_j$.

For the global system our unknown is
$
vv(U)=[vv(u)^((1)), vv(u)^((2)), ..., vv(u)^((K))]^"T"
$
// PC: is this right?!?!?
and we then gather our complete system
$
cal(M) dv(vv(U), t) = cal(L) vv(U) - tilde(vv(f)).
$
As we want to use standard ODE solvers, we want this on the standard ODE form, that is
$
dv(vv(U), t) = cal(M)^(-1)(cal(L) vv(U) - tilde(vv(f))).
$

= Numerical Experiments: Rotating Hill on a Square or Circular Mesh

== Testing the Solver

== Energy conservation
Energy is given by 
$
E(t) = 1/2 ||u(dot, t)||_(L^2(Omega))^2
$
which can be estimated by
$
E(t) approx 1/2 sum_(k=1)^K  integral_(Omega_k) lr(|u_h^((k)) (t, vv(x))|)^2 dif vv(x)
$
where
$
integral_(Omega_k) lr(|u_h^((k)) (t, vv(x))|)^2 dif vv(x)&=integral_(Omega_k) (sum_(n=1)^(N_p) u_n^((k))(t) phi.alt_n (vv(x)))^2 dif vv(x)\
&=integral_(Omega_k) sum_(i=1)^(N_p) sum_(j=1)^(N_p) u_i^((k))(t) u_j^((k))(t) phi.alt_i (vv(x)) phi.alt_j (vv(x)) dif vv(x)\
&= sum_(i=1)^(N_p) u_i^((k))(t) sum_(j=1)^(N_p) u_j^((k))(t) integral_(Omega_k) phi.alt_i (vv(x)) phi.alt_j (vv(x)) dif vv(x)\
&= sum_(i=1)^(N_p) u_i^((k))(t) [cal(M)vv(u)]_i \
&= (vv(u)^((k)))^"T" cal(M)^((k)) vv(u)^((k))
$
Using the global mass matrix, we get
$
 E(t) approx 1/2vv(u)^"T" cal(M) vv(u)
$

== Convergence Tests

= Appendix

We consider the toy problem $Omega = [0, 1]^2$
$
  nabla^2 u(x, y) = f(x, y), quad u(x, y) = g(x, y)
$
With solution
$
  u(x, y) = "bla"
$
Weak formulation
$
  integral_Omega nabla^2 u(x, y) phi(x, y) d x d y = integral_Omega f(x, y) phi(x, y) d x d y
$
Split into subdomains
$
  sum_(k=1)^K integral_(D^k) nabla^2 u(x, y) phi(x, y) d x d y = sum_(k=1)^K integral_(D^k) f(x, y) phi(x, y) d x d y
$
Greens identity
$
  integral_(D^k) nabla^2 u(x, y) phi(x, y) d x d y = - integral_(D^k) nabla u(x, y) nabla phi(x, y) d x d y + "boundary terms"
$
Consider
$
  Psi(mat(r ; s)) = mat(x ; y), quad J(r,s) = mat(x_r, x_s; y_r, y_s) = mat(pdv(x,r), pdv(x,s); pdv(y,r), pdv(y,s)), quad J(r, s) = J
$
Note that
$
  mat(x_r, x_s; y_r, y_s) mat(r_x, r_y; s_x, s_y) = mat(1, 0; 0, 1) arrow.r.double mat(x_r, x_s; y_r, y_s)^(-1) = mat(r_x, r_y; s_x, s_y)
$
Consider - this is not correct ? where is the transpose here !!!!
$
  nabla_((x,y)) u(x, y) &= mat(pdv(u,x) ; pdv(u, y)) = mat(pdv(r,x), pdv(r,y); pdv(s,x), pdv(s,y)) mat(pdv(u,r) ; pdv(u, s)) = mat(r_x, r_y; s_x, s_y) mat(pdv(u,r) ; pdv(u, s)) \ &= (mat(x_r, x_s; y_r, y_s)^(-1)) mat(pdv(u,r) ; pdv(u, s)) = (J^(-T)) nabla_((r,s)) u(r, s)
$
Consider - that is true
$
  integral_(D^k)u(x, y) d x d y = integral_I u(r, s) det(J(r, s)) d r d s
$
Consider - is it true?
$
  integral_(D^k) nabla_((x,y)) u(x, y) d x d y = integral_I (J^(-T)) nabla_((r,s)) u(r, s) det(J) d r d s
$
Mapping to canonical domain - is it true? (the functions are different, we should add tilde or something)
$
  I^"LHS"_k & = integral_(D^k) nabla_((x,y)) u(x, y) nabla_((x,y)) phi(x, y) d x d y \
   &= integral_I (J_k^(-T) nabla_((r,s)) u(r, s))^T (J_k^(-T) nabla_((r,s)) phi(r, s)) det(J_k) d r d s \
  &= integral_I nabla_((r,s)) u(r, s)^T (J_k^(-1) J_k^(-T)) nabla_((r,s)) phi(r, s) det(J_k) d r d s \
  &= integral_I nabla_((r,s)) u(r, s)^T (J_k^T J_k)^(-1) nabla_((r,s)) phi(r, s) det(J_k) d r d s
$
Consider
$
  I^"RHS"_k &= integral_(D^k) f(x, y) phi(x, y) d x d y = integral_I f(r, s) phi(r, s) det(J(r, s)) d r d s
$
Consider
$
  f(x, y) approx f_N (x, y) = sum_(i=1)^(N) hat(f_i) phi_i (x, y) = sum_(i=1)^(N) f_i h_i (x, y)
$
Based on Allan slides it should be the case that
$
  I^"LHS"_k &approx "diag"(J_k)[D_x^T M D_x + D_y^T M D_y] u_k \
  I^"RHS"_k &approx "diag"(J_k) M s_k
$
Note
$
  D_r = V_r V^(-1), quad D_s = V_s V^(-1)
$
Further
$
  D_x &= "diag"(bold(r)_x) D_r + "diag"(bold(s)_x) D_s \
  D_y &= "diag"(bold(r)_y) D_r + "diag"(bold(s)_y) D_s
$
The system of equations we are solving is the following
$
  - (J_k M)^(-1)[J_k D_x^T M D_x + J_k D_y^T M D_y] u_k = s_k
$

== Figuring out operations

Matrices proposed by chat
$ M_("sys") = 
mat(
  m_(11), m_(12), m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), m_(32), m_(33), m_(34) ;
  m_(41), m_(42), m_(43), m_(44) ;
) quad 
L_("sys") = mat(
  a_(11), a_(12), a_(13), a_(14) ;
  0, 0, 0, 0 ;
  a_(31), a_(32), a_(33), a_(34) ;
  a_(41), a_(42), a_(43), a_(44) ;
)
$
imposing BC
$
  u = mat(u_1, u_2, u_3, u_4) -> u_("new") = mat(u_1, g_2, u_3, u_4)
$
Consider
$
  L_("sys") u_("new") = mat(
  a_(11), a_(12), a_(13), a_(14) ;
  0, 0, 0, 0 ;
  a_(31), a_(32), a_(33), a_(34) ;
  a_(41), a_(42), a_(43), a_(44) ;
) mat(u_1 ;g_2; u_3; u_4) =\ mat(
  a_(11), 0, a_(13), a_(14) ;
  0, 0, 0, 0 ;
  a_(31), 0, a_(33), a_(34) ;
  a_(41), 0, a_(43), a_(44) ;
) mat(u_1 ;g_2; u_3; u_4) + mat(
  0, a_(12), 0, 0 ;
  0, 0, 0, 0 ;
  0, a_(32), 0, 0 ;
  0, a_(42), 0, 0 ;
) mat(u_1 ;g_2; u_3; u_4) = mat(u_1 ;g_2; u_3; u_4) = mat(
  a_(11), 0, a_(13), a_(14) ;
  0, 0, 0, 0 ;
  a_(31), 0, a_(33), a_(34) ;
  a_(41), 0, a_(43), a_(44) ;
) mat(u_1 ;g_2; u_3; u_4) +  mat(a_(12) g_2 ; bold(0); a_(32) g_2; a_(42) g_2)
$
Consider
$
  M_("sys") d/(d t) u_("new") = mat(
  m_(11), m_(12), m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), m_(32), m_(33), m_(34) ;
  m_(41), m_(42), m_(43), m_(44) ;
) mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4) =\ mat(
  m_(11), 0, m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), 0, m_(33), m_(34) ;
  m_(41), 0, m_(43), m_(44) ;
) mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4) + mat(
  0, m_(12), 0, 0 ;
  0, 0, 0, 0 ;
  0, m_(32), 0, 0 ;
  0, m_(42), 0, 0 ;
) mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4) = mat(
  m_(11), 0, m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), 0, m_(33), m_(34) ;
  m_(41), 0, m_(43), m_(44) ;
) mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4) + mat(m_(12)  d/(d t) u_2 ; bold(0); m_(32)  d/(d t) u_2 ; m_(42)  d/(d t) u_2
)
$
We have the system of equations with RHS:
$
  mat(
  m_(11), 0, m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), 0, m_(33), m_(34) ;
  m_(41), 0, m_(43), m_(44) ;
  ) mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4) + mat(m_(12)  d/(d t) u_2 ; bold(0); m_(32)  d/(d t) u_2 ; m_(42)  d/(d t) u_2) = mat("rhs"_1; d / (d t) g_2; "rhs"_3; "rhs"_4) -> d / (d t) u_2 = d / (d t) g_2
$

Consider
$
  f = mat(f_1, f_2, f_3, f_4), quad b = tilde(M)_("sys") f  = mat(
  m_(11), m_(12), m_(13), m_(14) ;
  0, bold(0), 0, 0 ;
  m_(31), m_(32), m_(33), m_(34) ;
  m_(41), m_(42), m_(43), m_(44) ;
) mat(f_1 ; f_2 ; f_3 ; f_4) = mat(b_1 ; 0 ; b_3 ; b_4)
$
The system becomes
$
  M_("sys") d/(d t) u_("new") = - L_("sys") u_("new") + tilde(M)_("sys") f + d / (d t) g_("bc")
$
Hence
$
mat(
  m_(11), 0, m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), 0, m_(33), m_(34) ;
  m_(41), 0, m_(43), m_(44) ;
)mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4)  + mat(m_(12) d/(d t) u_2  ;0 ; m_(32) d/(d t) u_2  ; m_(42) d/(d t) u_2)   =
 - mat(
  a_(11), 0, a_(13), a_(14) ;
  0, 0, 0, 0 ;
  a_(31), 0, a_(33), a_(34) ;
  a_(41), 0, a_(43), a_(44) ;
) mat(u_1 ;g_2; u_3; u_4) -  mat(a_(12) g_2 ; bold(0); a_(32) g_2; a_(42) g_2)
 + mat(b_1 ; 0 ; b_3 ; b_4) + mat(0; d / (d t) g_2; 0; 0)
$
We obtain the system as 5.2 in FEM book:
$
  mat(
  m_(11), 0, m_(13), m_(14) ;
  0, 1, 0, 0 ;
  m_(31), 0, m_(33), m_(34) ;
  m_(41), 0, m_(43), m_(44) ;
) mat(d/(d t) u_1 ; d/(d t) u_2 ; d/(d t) u_3 ; d/(d t) u_4) = - mat(
  a_(11), 0, a_(13), a_(14) ;
  0, 0, 0, 0 ;
  a_(31), 0, a_(33), a_(34) ;
  a_(41), 0, a_(43), a_(44) ;
) mat(u_1 ;g_2; u_3; u_4) + mat(b_1 ; 0 ; b_3 ; b_4) - mat(a_(12) g_2 ; bold(0); a_(32) g_2; a_(42) g_2) + mat( - m_(12)  d/(d t) g_2 ; d/(d t) g_2; - m_(32)  d/(d t) g_2 ; - m_(42)  d/(d t) g_2 )
$


#bibliography("references.bib")
