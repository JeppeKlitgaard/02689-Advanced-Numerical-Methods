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

We seek to represent the discrete trigonometric polynomial in terms of a nodal expansion over the interval $[0, 2π)$ using $N=2k : k ∈ ℕ$ points, $x_j$ for $j = 0, 1, ..., N-1$. Note the open end of the interval indicating that the last point, $x_(N-1)$ is not placed at the boundary but rather one grid spacing $h=(2π)/N$ inside the domain.

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

As a consequence of the periodic nature of the Fourier expansion, we find that any reconstruction in such a basis will have the property $u(0) = u(2π)$, allowing 
---

Derivation of first-order Fourier Differentatioin matrix D. \
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
NOOOO, what a mess
$
  h'_j (x) = 1/N ( N/2 cos(N/2 (x - x_j)) cot(1/2(x - x_j)) - 1/N sin(N/2 (x - x_j)) / (sin(1/2(x - x_j))^2))
$
consider $x_j = (2pi) / N j and x_k - x_j = (2pi) / N (k - j)$
$
  h'_j (x_k) = 1/N ( N/2 cos(2 pi (k - j)) cot(pi / N (k - j)) - 1/N sin(2 pi (k - j))) / (sin(pi / N (k - j))^2))
$



= Polynomial Methods <sec:polynomial>