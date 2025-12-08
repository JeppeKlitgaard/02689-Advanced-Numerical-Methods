#import "poster-template.typ": template
#import "@preview/dice:1.0.0"
#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval, TT, pdv
#import "preamble.typ": mref, mathformatter

#let rng_seed = dice.randominit()
#let authors = (
  (name: "Jeppe Klitgaard", email: "s250250@dtu.dk"),
  (name: "Tymoteusz Barcinski", email: "s221937@dtu.dk"),
  (name: "Pernille Christie", email: "s204249@dtu.dk"),
)
#let authors = dice.shuffle(authors, seed: rng_seed).at(0)
#set document(title: "Solving the Rotating Hill with Spectral Element Methods")
#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#show: template.with(authors)
#set math.equation(numbering: "(1)")

#set block(above: 0.8em, below: 0.8em)

// Column width: 240mm

// PC: Why so much whitespace above the line at the bottom of the document

= Rotating Hill Problem
We consider the conservative advection equation
$
pdv(u, t) + nabla dot (vv(a)u)=&0, quad (x,y)in Omega,thin t>0\
u(x,y,0) =& u_0(x,y) quad (x,y)in Omega
$ <eq:standard_eq>
where $vv(a)$ is some velocity field. We assume that the flow-field is divergence free, i.e.
$nabla dot vv(a) = 0.$
/* PC: exact solution might be a bit too much text
An exact solution to this problem is
$
u_"exact"(x,y,t) = exp(-((x-x_c (t))^2+(y-y_c (t)))/(2 sigma^2))
$*/

= Skew-Symmetric Split Form Advection
We use the chain rule to rewrite @eq:standard_eq, and average the two equations, which yields
$
pdv(u,t)+1/2[nabla dot (vv(a) u)+vv(a) dot nabla u]=0.
$ <eq:skew_symmetric>
We show below that this results in a skew-symmetric operator.

= Weak Formulations 
We divide the domain into smaller triangles, $Omega_k$, such that
$
Omega approx tilde(Omega) = union.big_(k=1)^K Omega_k.
$
Multiplying @eq:skew_symmetric by a test function $v$, integrating and using integration by parts, we get a weak form for the entire domain, where the boundary term disappears. We then require the residual to vanish in each domain, i.e.
#[
  #set math.equation(numbering: none)
  $
  integral_(Omega_k)v pdv(u,t) dif Omega_k + 1/2[&integral_(Omega_k) v vv(a) nabla u dif Omega_k &- integral_(Omega_k) u vv(a) nabla v dif Omega_k]=0.
  $
]
For elements along the boundary, we also get a term of 
$
1/2 integral_(partial Omega_j inter partial tilde(Omega)) v u(vv(a) dot vv(n) )dif s
$


= Upwind Numerical Flux

= Semi-Discrete Formulation
We divide $u$ into several elements $u^((k))$ such that $u(t, vv(x))=plus.o.big_k u^((k))(t, vv(x))$.
For each element, we approximate the function as 
$
u^((k))(t,vv(x)) approx u_h^((k))(t, vv(x))=sum_(n=1)^(N_p) u^((k))_n (t) phi.alt_n (vv(x)). quad
$
We are using a Galerkin method, meaning that $v=phi.alt_m$. We then combine this to a system for each element,
$
cal(M)^((k)) & dv(vv(u)^((k)),t)=cal(L)^((k))vv(u)^((k))
$
with
#[
  #set math.equation(numbering: none)
  $
  cal(L)^((k)) = -1/2 [ (A_x^((k)) - (A_x^((k)))^"T")  - (A_y^((k)) - (A_y^((k)))^"T") ]\
  A_x^((k)) = "diag"(a_x)cal(M)^((k))cal(D)_x, quad A_y^((k)) = "diag"(a_y) cal(M)^((k))cal(D)_y.
  $
]
/*A similar derivation also gives us the non-skew-symmetric system
$
cal(M)^((k)) & dv(vv(u)^((k)),t)= [& A_x  - A_y ]vv(u)^((k)).
$*/ // PC: We don't really use this - but should we?
A global system matrix is obtained by appropriately assembling the local matrices.

= Energy Stable Scheme // PC: is seperate section necessary?

Energy and its approximation is given by 
$
E(t) = 1/2 ||u(dot, t)||_(L^2(Omega))^2 approx 1/2(vv(u)(t))^"T" cal(M) vv(u)(t).
$
//PC: should we write u(t) here?
The scheme is energy-stable if $(dif E)/(dif t)< 0$, energy conserving if $(dif E)/(dif t)=0$ and energy unstable if $(dif E)/(dif t)>0$.

= Convergence Experiments

= Concluding Remarks ???

= References