#import "@preview/physica:0.9.4": curl, grad, tensor, pdv, dv, eval
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "@preview/cheq:0.2.2": checklist
#import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1"
#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#import "preamble.typ": mref, mathformatter

#show: checklist
#show: codly-init.with()

#codly(languages: codly-languages)

#set math.equation(numbering: "(1)")

#set heading(numbering: "1.a.1) ")
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

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: false, upright: true)
#set math.mat(delim:"[")

#let wider = h(3em)

#align(center,
  [#text(
  18pt,
  weight: "bold",
  [Assignment 1])
  #v(-6pt)
  #text(
    14pt,
    [02689 — Advanced Numerical Methods for Differential Equations
  ]
  )
])

Group XX:
- Jeppe Klitgaard <`s250250@dtu.dk`>
- Tymoteusz Barcinski <`s221937@dtu.dk`>
- Pernille Christie <`s204249@dtu.dk`>

= Fourier Spectral Methods <sec:fourier>

We are given the objective function
$
u'(x) ≔ 1/(2-cos(π x)) wider x ∈ [0, 2]
$

Which we manipulate to match the conventional domain $[0, 2π]$ of the complex Fourier coefficients $c_n$ in order to obtain 
$
u(x) ≔ 1/(2-cos(x)) wider x ∈ [0, 2π]
$ <eq:a_obj>

Where we recall the complex Fourier coefficients as:
$
c_n = 1/(2π) ∫_0^(2π) u(x) e^(-i n x) dif x wider n ∈ ℤ
$

Which for @eq:a_obj may be solved as:
$
c_n
&= 1/(2π) ∫_0^(2π) 1/(2-cos(x)) e^(-i n x) dif x wider n ∈ ℤ\
&= ...\
&= 1/(sqrt(3)(2 + sqrt(3))^(|n|))
$ <eq:a_coeffs>


Calculations:
$
c_(n+1) = 1/(2pi) integral_0^(2pi) cos((n+1)x) / (2 - cos(x)) dif x = 2 / (2pi) integral_0^(2pi) (cos(x) cos(n x)) / (2 - cos(x)) dif x - 1/(2pi)  integral_0^(2pi) cos((n-1)x) / (2 - cos(x)) dif x
$
We used the trigonometric identity
$
cos(alpha + beta) = cos(alpha) sin(beta) - sin(alpha) sin(beta) \
cos(alpha - beta) = cos(alpha) sin(beta) + sin(alpha) sin(beta) \
cos(alpha + beta) + cos(alpha - beta) = 2 cos(alpha) sin(beta)
$
Therefore
$
cos(n x + x) + cos(n x - x) = 2 cos(x) cos(n x)
$
Consider
$
integral_0^(2pi) (cos(x) cos(n x)) / (2 - cos(x)) dif x
&= integral_0^(2pi) ((2 - (2 -cos(x))) cos(n x)) / (2 - cos(x)) dif x\
&= 2 integral_0^(2pi) cos(n x) / (2 - cos(x)) dif x - integral_0^(2pi) cos(n x) dif x 
$
Note
$
2 c_n = 2 / (2 pi)  integral_0^(2pi) cos(n x) / (2 - cos(x)) dif x and integral_0^(2pi) cos(n x) dif x = 0 and c_(n-1) = 1/(2pi)  integral_0^(2pi) cos((n-1)x) / (2 - cos(x)) dif x
$
Therefore we obtain the characteristic equation
$
c_(n+1) = 4 c_n - c_(n-1)
$
If we just want to use induction we can now do
$
c_(n+1) = 4 / (sqrt(3)(2 + sqrt(3))^n) - (2 + sqrt(3)) / (sqrt(3)(2 + sqrt(3))^n) = (2 - sqrt(3)) / (sqrt(3)(2 + sqrt(3))^n) = 1 / (sqrt(3)(2 + sqrt(3))^(n+1))
$

If we go for the full solution without assuming that we know $c_n$. \
The characterisic polynomial is
$
alpha^2 - 4 alpha + 1 = 0 => alpha_1 = (4 + sqrt(12))/2 = 2 + sqrt(3) and alpha_2 = (4 - sqrt(12))/2 = 2 - sqrt(3)
$
Hence the solution is of the form
$
c_n = A_1 alpha_1 ^ n + A_2 alpha_2 ^n = A_1 (2 + sqrt(3))^n + A_2 (2 - sqrt(3))^n 
$
We need to establish initial condition $c_0, c_1$. From standard integraion rules
$
c_0 = 1/(2pi) integral_0^(2pi) cos(0) / (2 - cos(x)) dif x = 1/sqrt(3) and c_1 = 1/(2pi) integral_0^(2pi) cos(x) / (2 - cos(x)) dif x = 1/(sqrt(3)(2 + sqrt(3)))
$
Hence we need to solve the system
$
c_0 = A_1 + A_2 and c_1 = A_1 (2 + sqrt(3)) + A_2 (2 - sqrt(3)) 
$
The system has the solution
$
A_1 = 0 and A_2 = 1/sqrt(3)
$
Therefore
$
c_n = 1 / sqrt(3) (2 - sqrt(3))^n  = 1 / (sqrt(3) (2 + sqrt(3))^n)
$

== Convergence Behaviour of Truncated Fourier Expansion

Restricting the number of coefficients, $c_n$


This can be solved in a few ways:
1. Residue Theorem
2. Manipulate integrand to match standard integral (which standard integral though?)

== Convergence Behaviour of Discrete Fourier Coefficients


== Nodal Expansion

We seek to represent the discrete trigonometric polynomial in terms of a nodal expansion over the interval $[0, 2π)$ using $N=2k : k ∈ ℕ$ points, $x_j=(2π)/N j$ for $j = 0, 1, ..., N-1$. Note the open end of the interval indicating that the last point, $x_(N-1)$ is not placed at the boundary but rather one grid spacing $h=(2π)/N$ inside the domain.

To aide the reader through the treacherous ambiguity of how the nodal points are placed within the domain, see @fig:nodal_expansion.

#figure(
  cetz.canvas({
    import cetz.draw: *
    // Your drawing code goes here
    set-style(
      grid: (
        stroke: gray + 0.2pt,
        step: 1,
      ),
      line: (
        stroke: gray
      )
    )
  
    line((-2, 0), (2, 0))
    // Cross lines
    for i in (-2, -1, 0, 1, 2) {
      line((i, -6pt), (i, 6pt))
    }
    content((-2, 0.5), [$0$])
    content((2, 0.5), [$2π$])
  
    for i in (-2, -1, 0, 1) {
      circle((i, 0), radius: 3pt, fill: black)
    }
  
    content((-2.8, 0.0), [$x_j$])
    content((-2.8, -0.7), [$j$])
    for (k, i) in (-2, -1, 0, 1).enumerate() {
      content((i, -0.7), [$#k$])
    }
  }),
  caption: [Nodal Expansion. Note that last point $x_3$ is not located at the boundary, as this point is already given by $x_0$ due to periodicity of the trigonometric polynomials.]
) <fig:nodal_expansion>

As a consequence of the periodic nature of the Fourier expansion, we find that any reconstruction in such a basis will have the property $u(0) = u(2π)$.

We recall the Fourier series:
$
f(x) = sum_(-∞)^∞ c_n e^(i n x) wider x ∈ [0, 2π]
$ <eq:fourier_series>

Where the complex coefficients $c_n$ are given by:
$
c_n = 1/(2π) ∫_0^(2π) f(x) e^(-i n x) dif x
$ <eq:fourier_coeff>

We seek to find $c_n$ by evaluating the integral in @eq:fourier_coeff by quadrature using the trapezoidal rule:
$
∫_0^(2π) f(x) e^(-i n x) ≈ h sum_(j=0)^(N-1) ((f(x_(j+1) e^(-i n x_(j+1))) + f(x_j e^(-i n x_j)))/2)
$

Noting that the element $j=N-1$ of the sum contains the term $f(x_N) e^(-i n x_N) = f(x_0) e^(-i n x_0)$ by the aforementioned periodicity (which also holds for the factor $e^(-i n x)$ by inspection) we may rewrite the sum as follows upon realising that each term enters exactly twice, cancelling out the factor of $1/2$:
$
∫_0^(2π) f(x) e^(-i n x) ≈ h sum_(j=0)^(N-1) f(x_j) e^(-i n x_j)
$

Which enables the evaluation of the discretised complex coefficients $tilde(c)_n$ as:
$
tilde(c)_n
&= 1/(2π) (2π)/N sum_(j=0)^(N-1) f(x_j) e^(-i n x_j)\
&= 1/N sum_(j=0)^(N-1) f(x_j) e^(-i n x_j)
$ <eq:discrete_fourier_coeff>

Such that an objective function $f(x)$ may be approximated by expansion into a truncated Fourier series using the discrete coefficients in @eq:discrete_fourier_coeff as:
$
f(x)
&≈ sum_(n=-N/2)^(N/2) α_n tilde(c)_n e^(i n x)
&≈ sum_(n=-N/2)^(N/2) ( 1/N sum_(j=0)^(N-1) f(x_j) e^(-i n x_j) ) e^(i n x)\
&= sum_(j=0)^(N-1) f(x_j) underbrace(1/N sum_(n=-N/2)^(N/2) e^(-i n x_j) e^(i n x), h_j (x) )
$

Where
$
α_n ≔ cases(
  1/2 &quad |n| = N/2\
  1 &quad |n| < N/2
)
$

Is introduced to handle the overcounting arising from the fact that $e^(i N/2 x) = e^(- i N/2 x) med ∀ x=x_i$.

$h_j (x)$ is understood to scale the contributions of each node $x_j$ and may be expressed as:
$
h_j (x)
&= α_n/N sum_(n=-N/2)^(N/2) e^(-i n x_j) e^(i n x) = α_n/N sum_(n=-N/2)^(N/2) e^(i n θ) & θ ≔ x-x_j\
&= 1/N [(e^(i N/2 θ) + e^(-i N/2 θ))/2 + sum_(n=-N/2+1)^(N/2-1) e^(i n θ)]\
&= 1/N [sum_(n=-N/2)^(N/2) e^(i n θ) - (e^(i N/2 θ) + e^(-i N/2 θ))/2]
&= 1/N [sum_(n=-N/2)^(N/2) e^(i n θ) - cos(N/2 θ)]
$

Where the sum may be rewritten to closed form using the identity for a finite geometric series, $sum_(k=0)^M = a (1-r^M)/(1-r)$ with initial value $a=e^(-i N/2 θ)$ and common ratio $r=e^(i θ)$ with $M=N+1$ elements:
$
sum_(n=-N/2)^(N/2) e^(i n θ)
&= e^(-i N/2 θ) ⋅ (1-e^(i (N+1) θ))/(1- e^(i θ))\
&= e^(-i N/2 θ) ⋅ (e^(i θ/2) (e^(-i θ/2)-e^(i (N+1/2) θ)))/(e^(i θ/2) (e^(-i θ/2)- e^(i θ/2)))\
&= (e^(-i (1/2+N/2)θ)-e^(i (1/2+N/2) θ))/(e^(-i θ/2)- e^(i θ/2))\
&= (2i ⋅ sin((N+1)/2 θ))/(2i ⋅ sin(θ/2))
= sin((N+1)/2 θ)/sin(θ/2)
$

Which in turn enables a final rewrite of $h_j (x)$ using the identity $sin(α ± β) = sin(α) cos(β) ± cos(α)sin(β)$:
$
h_j (x)
&= 1/N [sin((N+1)/2 θ)/sin(θ/2) - cos(N/2 θ)]\
&= 1/N [(sin(N/2 θ) cos(θ/2) + cos(N/2 θ) sin(θ/2))/sin(θ/2) - cos(N/2 θ)]\
&= 1/N [sin(N/2 θ) cot(θ/2) + cos(N/2 θ) - cos(N/2 θ)]\
&= 1/N sin(N/2 (x-x_j)) cot((x - x_j)/2) &qed
$ <eq:lagrange_polynomials>

Inspecting @eq:lagrange_polynomials easily reveals that the polynomials take the form of the Dirac delta function over the nodes, $h_j (x_i) = δ_(i j) quad ∀ i, j ∈ {0, ..., N-1}$.



== Derivation of first-order Fourier Differentatioin matrix $D$
$
  I_N u(x) = sum_(j=0)^N u(x_j) h_j (x) and d/(d x) I_N u(x) = sum_(j=0)^N u(x_j) h'_j (x)
$
where
$
  h_j (x) = 1/N sin(N/2 (x - x_j)) cot(1/2(x - x_j)) = 1/N (sin(N/2 (x - x_j)) cos(1/2(x - x_j))) / (sin(1/2(x - x_j)))
$
consider
$
d / (d x) sin(N/2 (x - x_j)) cos(1/2(x - x_j)) = \
N/2 cos(N/2 (x - x_j)) cos(1/2 (x - x_j)) - sin(N/2 (x - x_j)) sin(1/2 (x - x_j)) \
d / (d x) sin(1/2(x - x_j)) = 1/2 cos(1/2(x - x_j))
$
hence
$
  h'_j (x) = (d / (d x) (sin(N/2 (x - x_j)) cos(1/2(x - x_j)))) / (sin(1/2(x - x_j))) - (sin(N/2 (x - x_j)) cos(1/2(x - x_j))^2) / (sin(1/2(x - x_j))^2)
$
NOOOO, what a mess \
This works
$
  h'_j (x) = 1/N ( N/2 cos(N/2 (x - x_j)) cot(1/2(x - x_j)) - 1/N sin(N/2 (x - x_j)) / (sin(1/2(x - x_j))^2))
$
consider $x_j = (2pi) / N j and x_k - x_j = (2pi) / N (k - j)$
$
  h'_j (x_k) = 1/N ( N/2 cos(pi (k - j)) cot(pi / N (k - j)) ) - 1/N sin(pi (k - j)) / (sin(pi / N (k - j))^2) = \
  1/2 (-1)^(k-j) cot(pi / N (k - j)) - 0 = 1/2 (-1)^(k-j) cos(pi / N (k - j)) / sin(pi / N (k - j))
$



= Polynomial Methods <sec:polynomial>