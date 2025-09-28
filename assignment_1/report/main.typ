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
  [Assignment 1])
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

= Fourier Spectral Methods <sec:fourier>

== Truncated Fourier Expansion <sec:1a>
We are given the objective function
$
tilde(u)(x) ≔ 1/(2-cos(π x)) wider x ∈ [0, 2]
$

Which we manipulate to match the conventional domain $[0, 2π]$ of the complex Fourier coefficients $c_n$ in order to obtain 
$
u(x) ≔ 1/(2-cos(x)) wider x ∈ [0, 2π]
$ <eq:a_obj>

=== Closed-Form Analytical Fourier Coefficients

In order to determine a closed-form solution to the Fourier coefficients of a Fourier series expansion of @eq:a_obj, we first recall the complex Fourier coefficients as:
$
c_n
&= 1/(2π) ∫_0^(2π) u(x) e^(-i n x) dif x wider n ∈ ℤ\
&= 1/(2π) ∫_0^(2π) 1/(2-cos(x)) e^(-i n x) dif x\
$ <eq:fourier_coeff_unsolved>

Inspection of @eq:a_obj immediately reveals that $u(x)$ is an _even_ and _real_ function, which in turn implies that the sine component of the integrand must vanish. This can be shown by using Euler's formula to express @eq:fourier_coeff_unsolved as:
$
c_n
&= 1/(2π) [∫_0^(2π) cos(n x)/(2-cos(x)) dif x + cancel(∫_0^(2π) (i sin(n x))/(2-cos(x)) dif x)]\
&= 1/(2π) ∫_0^(2π) cos(n x)/(2-cos(x)) dif x
$ <eq:fourier_coeff_unsolved_2>

Where the cancellation becomes obvious upon realising that the integrand is odd centered around $π$ (as it is a quotient of an even and an odd function) and thus integrates to zero over the interval $[0, 2π]$.

We seek to find a closed-form, analytical solution of @eq:fourier_coeff_unsolved_2 by means of induction by considering:
$
c_(n+1) = 1/(2pi) integral_0^(2pi) cos((n+1)x) / (2 - cos(x)) dif x = 2 / (2pi) integral_0^(2pi) (cos(x) cos(n x)) / (2 - cos(x)) dif x - 1/(2pi)  integral_0^(2pi) cos((n-1)x) / (2 - cos(x)) dif x
$ <eq:1a_inductive_step_1>

Where the latter expression arises from the use of the trigonometric identities
$
cos(alpha ± beta) = cos(alpha) cos(beta) ∓ sin(alpha) sin(beta)\
⇕\
cos(n x + x) + cos(n x - x) = 2 cos(x) cos(n x)
$

We consider the following integral:
$
∫_0^(2π) (cos(x) cos(n x)) / (2 - cos(x))
&= ∫_0^(2π) ((2 - (2 -cos(x))) cos(n x)) / (2 - cos(x)) dif x\
&= 2 ∫_0^(2π) cos(n x) / (2 - cos(x)) dif x - ∫_0^(2π) cos(n x) dif x\
&= 4π  c_n - 2π δ_(0 n)
$

Such that the inductive step from @eq:1a_inductive_step_1 becomes:
$
c_(n+1) = 4c_n - 2δ_(0n) - c_(n-1)
$ <eq:1a_inductive_step_2>

We consider the case $n>0$ such that @eq:1a_inductive_step_2 becomes:
$
c_(n+1) = 4c_n - c_(n-1)
$

Yielding the characteristic polynomial:
$
α^2-4α+1=0
$

With solutions
$
α_(1,2) = 2 ± sqrt(3)
$

Leading to coefficients $c_(n)$ of the form:
$
c_(n+1) = A_1 α_1^n + A_2 α_2^n = A_1 (2 + sqrt(3))^n + A_2 (2 - sqrt(3))^n
$

Where $A_1, A_2$ may be determined by considering cases $c_0, c_1$:
$
c_0 &= A_1 + A_2 &&= 1/(2π) ∫_0^(2π) cos(0)/(2-cos(x)) dif x &&= 1/sqrt(3)\
c_1 &= A_1(2+sqrt(3)) + A_2(2-sqrt(3)) &&=1/(2π) ∫_0^(2π) cos(x)/(2-cos(x)) dif x &&= 1/(sqrt(3)(2+sqrt(3)))\
$ <eq:1a_initial_steps>
#math.equation(block: true, numbering: none)[
	$
⇕
$	
]
$
A_1 = 0 wider and wider A_2 = 1/sqrt(3)
$

Thus providing a closed-form solution for the coefficients $c_n$ for the case $n>0$:
$
c_n = 1/sqrt(3) (2-sqrt(3))^n
$

We may generalise this solution to the $n=0$ case by simply observing that $c_0 = 1/sqrt(3)(2-sqrt(3))^0 = 1/sqrt(3)$ as in @eq:1a_initial_steps.

The case $n<0$ may be elucidated by symmetry arguments using the fact that the function $u(x)=1/(2-cos(x))$ is even, which implies that the Fourier coefficients themselves will be even: $c_n = c_(-n)$.
This allows the complex Fourier coefficients $c_n$ to be given as a single closed-form solution:
$
c_n = 1/sqrt(3)(2-sqrt(3))^(|n|) = 1/(sqrt(3) (2+sqrt(3))^(|n|)) wider n ∈ ℤ
$

We thus find the asymptotic decay rate of the continuous Fourier coefficients to be exponential and of order $cal(O)((2+sqrt(3))^(-|n|))$.


=== Convergence Behaviour of Truncated Fourier Expansion

We may then investigate the truncation error $||u(x) -cal(P)_N u(x)||_2$ where $cal(P)_N$ denotes a truncated Fourier expansion with $N$ modes.

We simply evaluate the exact function and the truncated Fourier expansion on a fine grid over $x∈[0, 2π]$ for $N∈{1, ..., 50}$ to obtain the graph shown in @fig:1a_fourier_truncation_convergence.

#figure(
  image("output/1a_fourier_truncation_convergence.png"),
  caption: [Truncation Error $||u(x) -cal(P)_N u(x)||_2$ for $u(x) = 1/(2-cos(x))$ over a range of truncation levels, $N$]
) <fig:1a_fourier_truncation_convergence>

We observe that we reach machine precision at $N ≈ 27$, after which the truncation error is not improved by employing more modes in the expansion.

== Convergence Behaviour of Discrete Fourier Coefficients <sec:1b>

After having found a closed-form expression for the analytical Fourier coefficients in @sec:1a, we would like to compare this with the discrete Fourier coefficients, $tilde(c)_n$, which may be obtained using a Fast Fourier Transform (FFT).

The discrete Fourier coefficients, $tilde(c)_n$, are derived below in @sec:1c1 as part of the even nodal expansion, and may be represented as:
$
tilde(c)_n = 1/N sum_(j=0)^(N-1) f(x_j) e^(-i n x_j)
$

These may be computed using the _Discrete Fourier Transform Method_ by leveraging the FFT routines found in modern scientific computing packages such as `numpy` and `scipy`. Care needs to be taken around the normalisation and index conventions, which is handled using `norm="forward"` which scales the transform by $1 \/ N$, and the appropriate `fftshift` operations.

We seek to investigate the effect of aliasing and the errors it introduces in the coefficients. Aliasing refers to the phenomenon where high-frequency modes in a continuous signal are indistinguishable from lower-frequency modes when sampled discretely. Effectively, the higher-frequency modes are 'folded' or 'aliased' into the lower-frequency components, and the resulting discrete coefficient $tilde(c)_n$ becomes the sum of the true coefficient $c_n$ and all of its higher-frequency aliases $c_(k±m N)$. We should be able to observe this effect as a deviation between the discrete and analytical coefficients, particularly at higher wavenumbers. Aliasing effects are often likened to wheels going in reverse in video recordings.

We compute the discrete Fourier coefficients with different levels of discretisation, $N∈{4, 8, 16, 32, 64}$ and compare them against the analytical coefficients $c_n$ for which we obtained an analytical expression in @sec:1a. As seen in @fig:1b_errors, we observe that the coefficients of $tilde(c)_n$ decay exponentially with respect to $|n|$ as predicted by theory. The effect of aliasing is greatest at the higher frequency modes, as $n → |N/2|$. While this may initially be concerning, we note that the coefficients $c_n$ decay faster for the higher-frequency modes for any sufficiently smooth function, as seen here for $u(x)$.

#figure(
  image("output/1b_fourier_dft_error_analysis.png"),
  caption: [Convergence and errors associated with aliasing effects in discretisation of Fourier coefficients for $u(x) = 1/(2-cos(x))$.]
) <fig:1b_errors>


From @fig:1b_errors it becomes clear that $N$ must be chosen sufficiently high to capture the significant modes of the underlying function, which matches our expectation that the smoothness of the function sets a lower bound on $N$ required to achieve a given tolerance. Any frequencies above the $N\/2$#super("th") mode, corresponding to the _Nyquist Frequency_ are not captured in the expansion and instead enter in the lower frequency modes as an error. As such, we can summarise the smoothness condition as: _sufficient smoothness is implied by the signal having no significant contributions above the Nyquist Frequency_.

In the lower panel of @fig:1b_errors we find spectral convergence for a chosen mode, $c_0$, reaching machine precision at a discretisation level of $N=32$.

As such, we find that the trigonometric polynomials form a good basis in which to discretise a signal for any sufficiently smooth, periodic function and can be used to faithfully represent such a signal provided it is without significant contributions above the _Nyquist Frequency_, which is determined by the chosen number of modes, $N$.

== Nodal Expansion <sec:1c>

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

=== Closed-Form Lagrange Polynomial Derivation <sec:1c1>

We seek to derive a closed-form solution for the Lagrange polynomials using a nodal expansion of the discrete trigonometric polynomial by first recalling the Fourier series:
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
&≈ sum_(n=-N/2)^(N/2) α_n tilde(c)_n e^(i n x)\
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
is introduced to handle the overcounting arising from the fact that $e^(i N/2 x) = e^(- i N/2 x) med ∀ x=x_i$.

$h_j (x)$ is understood to scale the contributions of each node $x_j$ and may be expressed as:
$
h_j (x)
&= α_n/N sum_(n=-N/2)^(N/2) e^(-i n x_j) e^(i n x) = α_n/N sum_(n=-N/2)^(N/2) e^(i n θ) & θ ≔ x-x_j\
&= 1/N [(e^(i N/2 θ) + e^(-i N/2 θ))/2 + sum_(n=-N/2+1)^(N/2-1) e^(i n θ)]\
&= 1/N [sum_(n=-N/2)^(N/2) e^(i n θ) - (e^(i N/2 θ) + e^(-i N/2 θ))/2]
&= 1/N [sum_(n=-N/2)^(N/2) e^(i n θ) - cos(N/2 θ)]
$

Where the sum may be rewritten to closed form using the identity for a finite geometric series, 
$
sum_(k=0)^M a r^k = a (1-r^M)/(1-r),
$ 
with initial value $a=e^(-i N/2 θ)$ and common ratio $r=e^(i θ)$ with $M=N+1$ elements:
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

=== Visualisation of Lagrange Polynomials

Having now obtained an expression for the Lagrange Polynomials, $h_j (x)$ as given in @eq:lagrange_polynomials, we visualise the polynomials for the $N=6$ case as shown in @fig:lagrange_polynomials.

#figure(
  image("output/1c_lagrange_polynomials.png"),
  caption: [Lagrange Polynomials $h_j (x)$ for $N=6$. Note that $h_j (x_i) = δ_(i j)$ where $x_i$ are demarked by the dashed vertical lines.]
) <fig:lagrange_polynomials>


=== Fourier Differentiation Matrix
To derive the first-order Fourier Differentiation matrix $D$, we differentiate the Fourier Interpolation function
$
  I_N u(x) = sum_(j=0)^(N-1) u(x_j) h_j (x) "where" h_j (x) = 1/N sin(N/2 (x - x_j)) cot(1/2(x - x_j))
$
We seek the following expansion
$
  d/(d x) I_N u(x) = sum_(j=0)^(N-1) u(x_j) h'_j (x)
$

// We therefore consider the derivative of $h_j$, using the product rule, 
// $
// d / (d x) sin(N/2 (x - x_j)) cos(1/2(x - x_j)) = \
// N/2 cos(N/2 (x - x_j)) cos(1/2 (x - x_j)) - sin(N/2 (x - x_j)) sin(1/2 (x - x_j)) \
// d / (d x) sin(1/2(x - x_j)) = 1/2 cos(1/2(x - x_j))
// $
// hence
// $
//   h'_j (x) = (d / (d x) (sin(N/2 (x - x_j)) cos(1/2(x - x_j)))) / (sin(1/2(x - x_j))) - (sin(N/2 (x - x_j)) cos(1/2(x - x_j))^2) / (sin(1/2(x - x_j))^2)
// $
// NOOOO, what a mess \

The derivative of Lagrange polynomial located at point $x_j = (2pi) / N j$ is
$
  h'_j (x) = 1/N ( N/2 cos(N/2 (x - x_j)) cot(1/2(x - x_j)) - 1/N sin(N/2 (x - x_j)) / (sin(1/2(x - x_j))^2))
$
We consider the equidistant grid of the same points that was used to construct Lagrange polynomials, namely $x_k = (2pi) / N k$, resulting in $x_k - x_j = (2pi) / N (k - j)$. Inserting into the above expression leads to the following simplification
$
  h'_j (x_k) = 1/N ( N/2 cos(pi (k - j)) cot(pi / N (k - j)) ) - 1/N sin(pi (k - j)) / (sin(pi / N (k - j))^2) = \
  1/2 (-1)^(k - j) cot(pi / N (k - j)) - 0 = 1/2 (-1)^(k - j) cos(pi / N (k - j)) / sin(pi / N (k - j)).
$
Note that $(-1)^(k - j) = (-1)^(k + j)$. Thus, the entries of the Fourier Differentiation matrix $D$ are
$
D_(k j) = h'_j (x_k) = cases(
  1/2(-1)^(k + j)cot(pi / N (k - j)) &"if" j != k, 
  0 &"if" j=k
)
$
// TODO: There is something off with our signs compared to the book.

 // skew-symmetry D_{nj} = −D_{jn},

== Fourier Differentiation Routine
// Comment on the fact that the function is $C^(infinity)$ therefore we get spectral convergence.

We consider the function
$
  v(x) = exp(sin(x)) quad x in [0, 2pi]
$
where we changed the argument to match the interval of the Fourier Differentiation matrix D. The derivative of $v$ is
$
  v'(x) = exp(sin(x)) cos(x) quad x in [0, 2pi]
$

#figure(
  image("output/1d_convergence_test.png"),
  caption: [Truncation error of $max_(j in [0, N-1]) {|v'(x_j) - D v(x_j)|}$ across different $N$ in the semi-log plot.]
) <fig:1d_convergence_test>

@fig:1d_convergence_test presents the convergence test for increasing $N$. We observe spectral convergence which is justified by the fact that $v in C^(infinity)([0, 2pi])$, which implies that function $v$ is smooth and hence Fourier Series converges rapidly.
// and uniformly.

// MORE MATH JUSTIFICATION NEEDED.

== Convergence in the $L^2$-norm

We consider the sequence of functions $w^i (x)$ with changed argument such that the domain of each function is $[-2pi, 2pi]$. Consider the analytical integration
$
  w^(i+1)(x) = integral_(-2pi)^(2pi) w^i (x) d x quad "such that" quad w^(i+1) in C^(i)([-2pi, 2pi]) quad "where" quad w^i (x) = (d w^(i+1)) / (d x) (x)
$
We obtain the following functions
$
  w^0 (x) = cases(
    -cos(x) &"if" x in [-2pi, 0), 
    cos(x) &"if" x in [0, 2pi]
  ) quad
  w^1 (x) = cases(
    -sin(x) &"if" x in [-2pi, 0), 
    sin(x) &"if" x in [0, 2pi]
  ) 
$
$
  w^2 (x) = cases(
    cos(x) - 1 &"if" x in [-2pi, 0), 
    -(cos(x) - 1) &"if" x in [0, 2pi]
  ) quad
  w^3 (x) = cases(
    sin(x) - x &"if" x in [-2pi, 0), 
    -(sin(x) - x) &"if" x in [0, 2pi]
  )
$
@fig:1e_functions_wi shows the functions $w^i (x) "for" x in [-2pi, 2pi]$. Note that the function $w^0$ is not periodic on the domain and has a discontinuity at $x=0$, whereas the function $w^i "for" i >= 1$ are periodic on $[-2pi, 2pi]$.

#figure(
  image("output/1e_functions.png", width: 80%),
  caption: [Functions $w^i (x) "for" x in [-2pi, 2pi]$]
) <fig:1e_functions_wi>

Consider the derivative of the interpolating Lagrange polynomial
$
  I_N (w^(i))' (x) = sum_(j=0)^N w^(i)(x_j) h'_j (x)
$
We can write the expression
$
w^(i-1)(x) = (w^(i))' (x) = sum_(j=0)^N w^(i)(x_j) h'_j (x) + tau_N (x)
$
Consider the $L^2$ norm of the truncation error
$
||tau_N||^2_(L^2) = integral_(-2 pi)^(2pi)(tau_N (x))^2 d x = integral_(-2 pi)^(2pi)(w^(i-1)(x) - sum_(j=0)^N w^(i)(x_j) h'_j (x))^2 d x
$
Consider the estimate of the $L^2$ norm by the discrete $L^2$ norm obtained via a quadrature rule with $M$ points
$
  ||tau_N||^2_(L^2) = integral_(-2 pi)^(2pi)(tau_N (x))^2 d x = sum_(m=0)^(M) (tau_N (x_m))^2 w_m + eta_M
$
where $eta_M$ is the truncation error of the quadrature rule, and we use the formula for quadrature rule for periodic functions. We use the quadrature routine from `scipy` with the number of points $M$ such that the truncation error $eta_M$ is below machine precision. This allows us to observe only the truncation error $tau_N$ for different $N$ when computing the discrete $L^2$ norm on the convergence test. Note that this is a generalization compared to the Fourier Differentiation matrix D, which constructs the Lagrange polynomials on a specific uniform grid and evaluates the derivative on the same uniform grid. In the above formulation, the interpolating Lagrange polynomial is constructed on the uniform grid with $N$ points, but the derivative is evaluated on the grid of size $M$.


@fig:1e_convergence_test presents the convergence test for discrete first derivative of $w^i$. We have the following
$
  w^1 in C^(0)([-2pi, 2pi]) quad w^2 in C^(1)([-2pi, 2pi]) quad w^3 in C^(2)([-2pi, 2pi])
$
Therefore the upper bound for the algebraic convergence is $1 "for" w^1$, $2 "for" w^2$, and $3 "for" w^3$. We observe that for derivatives of $w^2, w^3$ the convergence rate is higher.
Note, that since the function $w^0 = (w^1)'$ does not have a continuous first derivative and is not periodic the interpolating Lagrange polynomials do not converge point wise for points $x in {-2pi, 0, 2pi}$ and the neighboring points, and the Gibbs oscillations can be observed. However the convergence in the $L^2$ norm is achieved although only linear based on the convergence test.

// COMMENT: $w^0$ should converge in the L2 space theoretically. 

#figure(
  image("output/1e_convergence_test_wi.png"),
  caption: [Convergence test in the log-log plot for approximating the derivative of $w^i$ with Lagrange polynomials in a loglog plot.]
) <fig:1e_convergence_test>


// TODO: CHANGE FIGURE FROM GIT
// #figure(
//   image("output/1f_FFT_D_CPU_time.png"),
//   caption: [TODO]
// ) <fig:1e_functions_wi>

// Consider $N = 5$
// $
// (d f) / (d x) (x) = sum_(j=0)^N f(x_j) (d h_j) / (d x) (x)
// $

// Consider
// $
//   f(x) = sum_(j=0)^N f_j h_j (x) = bold(f)^T bold(h) (x) + tau(x)
// $
// Further
// $
//   ||f||^2_N = integral (bold(f)^T bold(h) (x) + tau(x))(bold(f)^T bold(h) (x) + tau(x)) d x = \
//   integral (bold(f)^T bold(h) (x))(bold(f)^T bold(h) (x)) d x + 2 integral bold(f)^T bold(h) (x) tau(x) d x + integral tau^2(x) d x = \
//   integral (bold(f)^T bold(h) (x))(bold(f)^T bold(h) (x)) d x + tilde(tau) = bold(f)^T (V V^T)^(-1) bold(f) + tilde(tau)
// $
// Consider
// $
//   ||tau||^2_N = ||f - bold(f)^T bold(h)||^2_N = integral (f - bold(f)^T bold(h))(f - bold(f)^T bold(h)) d x = \
//   integral (f(x))^2 d x - 2 integral (f(x) bold(f)^T bold(h) (x)) d x + integral (bold(f)^T bold(h) (x))(bold(f)^T bold(h) (x)) d x = \
//   ||f||^2_N - 2 integral (f(x) bold(f)^T bold(h) (x)) d x + bold(f)^T (V V^T)^(-1) bold(f)
// $
// Alan claim blackboard
// $
//   ||e||^2_N = sqrt(e^T M e)
// $
// Again
// $
//   ||tau||^2_N = integral (bold(tau)^T bold(h) (x) + eta(x))(bold(tau)^T bold(h) (x) + eta(x)) d x = \
//   ||eta||^2_N - 2 integral (eta(x) bold(tau)^T bold(h) (x)) d x + bold(tau)^T (V V^T)^(-1) bold(tau)
// $



== Numerical Differentiation using Fast Fourier Transform

The derivative of function $v$ can be computed via FFT by the following formula
$
  d / (d x) P_N v(x) = sum_(n=0)^(N-1) c_n (i n) e^(i n x)
$
Therefore we proceed by obtaining the coefficients $tilde(c)_n$ by the FFT across the equidistant grid of $N$ points, and subsequently we compute the inverse FFT of the set of coefficients ${i n c_n}_(n=0)^(N-1)$. The comparison of CPU time between the computation of derivative of function $v$ via the Fourier Differentiation matrix D and the FFT are presented in @fig:1f_cpu. The break point beyond which the FFT is faster than Fourier Differentiation matrix D was found to be $N=10^3$ on MacBook M4 Pro. The asymptotic convergence rate for Fourier Differentiation matrix D was found to be $O(N^2)$ and the asymptotic convergence rate for FFT was found to be $O(N log(N))$ for $N > 10^2$ which is consistent with the theory.

Note that for the comparison the function $v in C^(infinity)([0, 2pi])$ was used and spectral convergence is observed on the right plot of @fig:1f_cpu - The level of machine precision is achieved with less than $N=50$. Therefore for this function the differentiation via the Fourier Differentiation matrix D should be preferred. Functions which have only up to specific number of continuous derivatives converge algebraically and would require higher $N$ to achieve machine precision.

Generally, the results imply that the choice between Fourier Differentiation matrix D and FFT for derivative computation should be based on the expected number of required $N$, while the number of required $N$ to achieve desired precision should be based on the smoothness of the function to be differentiated.

#figure(
  image("output/1f_FFT_D_CPU_time.png"),
  caption: [Comparison of CPU time of Fourier Differentiation matrix D and the FFT for derivative computation of function $v$ for increasing $N$ in a loglog plot.]
) <fig:1f_cpu>


= Polynomial Methods <sec:polynomial>
#counter(heading).update((2, 7))

== Jacobi Polynomials
The Jacobi polynomials are given as a recursion,
$
P_0^((alpha,beta))(x)&=1\
P_1^((alpha,beta))(x)&=1/2(alpha-beta+(alpha+beta+2)x)\
P_n^((alpha,beta))(x)&=1/(a_(n,n-1)^((alpha,beta)))((a_(n-1,n-1)^((alpha, beta))+x)P_(n-1)^((alpha,beta))-a_(n-2,n-1)^((alpha, beta))P_(n-2)^((alpha,beta))),
$ <eq:jacobi_polynomials>
where
$
a_(n-1,n)^((alpha,beta))&=(2(n+alpha)(n+beta))/((2n+alpha+beta+1)(2n+alpha+beta))\
a_(n,n)^((alpha,beta))&=(alpha^2-beta^2)/((2n+alpha+beta+2)(2n+alpha+beta))\
a_(n+1,n)^((alpha,beta))&=(2(n+1)(n+alpha+beta+1))/((2n+alpha+beta+2)(2n+alpha+beta+1))
$
for $n>=0$, and $a^((alpha,beta))_(-1,0)=0$. 

To calculate the Jacobi polynomials, we have to calculate the recursion with the correct $a$-values. The code, `jacobi_p`, is implemented in such a way that the polynomials are calculated "from the bottom up", to prevent the costs of recursion. To check our implementation, we plot the first six Lagrange and Chebyshev polynomials, see @fig:jacobi_polynomials. For the Chebyshev polynomials, we note that they are scaled Jacobi polynomials,

$
T_n^((-1/2, -1/2)) (x)=(Gamma(n+1)Gamma(1/2))/Gamma(n+1/2) P_n^((-1/2, -1/2)) (x).

$
The plots align with the plots given in lectures, so we assume that our correction is correct.

#figure(
  image("output/2h_jacobi_polynomials.png"),
  caption: [Chebyshev polynomials $T_n^((-1/2, -1/2)) (x)$ and Legendre polynomials $P_n^((0,0)) (x)$ for $n in {0,...,5}$ calculated by the `jacobi_p` routine. Note that for the Chebyshev polynomial $P_n^((-1/2, -1/2)) (x)$ has been scaled.]
) <fig:jacobi_polynomials>

== Numerical Experiments
We calculate the abscissas $x_j$ using translated versions of the functions JacobiQP and JacobiGL, given on Learn. We then calculate the first $K=200$ coefficients as follows, 
$
  tilde(u)_k = 1/gamma_k sum_(j=0)^N u(x_j) phi_k (x_j) w_j, wide "for" k = 0,...,K,
$ <eq:nodal_values>
where in our case $u(x) = 1/(2 - cos(pi (x+1))$, $x in  [-1,1]$ and we use the Legendre polynomials, i.e. $phi_k (x) = P^(0,0)_k (x)$. Notice that we have redefined $u$ such that it covers the correct interval. A plot of the calculated coefficients can be seen in @fig:calculated_coefficients. 

#figure(
  image("output/2i_coefficients.png", width: 80%),
  caption: [Plot of calculated coefficients $|c_k|=|tilde(f)_k|$ for $k = 0,...,K$ and varying N in a semi-log plot. ]
) <fig:calculated_coefficients>

We see that for $N=K$ the coefficients behave as expected, where they decay towards 0 as $k$ grows. For $N < K$, we see that the coefficients follows the expected pattern until the error grows again. We calculate the interpolant as 
$
I_N u(x)=sum_(k=0)^(K) hat(u)_k thin phi_k(x),
$
which allows us to examine the interpolation error. Our plot of the error, see @fig:error_DFT, shows the error we expect from the coefficients, namely that we only get satisfactory results when $N=K$.

#figure(
  image("output/2i_error_DTM.png", width: 80%),
  caption: [Truncation error $||u(x) -I_N u(x)||_2$ of the Discrete Polynomial Transform for the function $tilde(u)(x) = 1/(2 - cos(pi (x+1))$ using 200 coefficients and for varying grid sizes in a semi-log plot.]
  
) <fig:error_DFT>

To explain this error, we look at how we calculate the coefficients,
$
tilde(f)_k=&(f,phi_k)_N/(||phi_k||_N^2)=(sum_(m=0)^infinity hat(f)_m phi_m,phi_k)_N/(||phi_k||_N^2) \
=&(sum_(m=0)^N hat(f)_m phi_m,phi_k)_N/(||phi_k||_N^2)+(sum_(m=N+1)^infinity hat(f)_m phi_m,phi_k)_N/(||phi_k||_N^2)\
=&1/(||phi_k||_N^2)sum_(m=0)^N hat(f)_m (phi_m,phi_k)_N+1/(||phi_k||_N^2) sum_(m=N+1)^infinity hat(f)_m (phi_m,phi_k)_N
$

For $k<=N$, we get
$
tilde(f)_k=hat(f)_k +1/(||phi_k||_N^2) sum_(m=N+1)^infinity hat(f)_m (phi_m,phi_k)_N,
$
due to the orthogonality of the $phi$-basis. This means that we have the exact coefficient plus some error. We know that $(phi_m,phi_k)_N in PP^(m+k)$ and that the quadrature is exact for any $f in PP^(2N-1)$, so we get errors from the sum term when $m+k>2N-1$. This is also true for $k>N$, but in this case, we no longer have the exact term and more terms will have a too high degree, resulting in large errors. This matches the image we see in @fig:calculated_coefficients.

== Generalized Vandermonde Matrix <sec:polynomial_j>
To construct the generalized Vandermonde matrix, we find polynomials and abscissas as in the previous exercise and gather in a matrix
$
bold(cal(V)) = mat(
  phi_1(x_1), phi_2(x_1), ..., phi_N (x_1);
  phi_1(x_2), phi_2(x_2), ..., phi_N (x_2);
  dots.v, dots.v, dots.down, dots.v;
  phi_1(x_N), phi_2(x_N), ..., phi_N (x_N)
)
$
To calculate the Lagrange polynomials on a uniform grid, $tilde(x)_i$ for $i=0,...,99$, we solve the equation 

$
bold(cal(V))^top bold(h)(tilde(x)_i)=bold(phi)(tilde(x)_j)
$
for each $i=0,...,99$ in the uniform grid. The resulting polynomials can be seen in @fig:lagrange_polynomials_2. We see that the Lagrange polynomials evaluated at the abscissas take the expected values of 1 or 0, specifically that $h_j (x_i) = δ_(i j)$, indicating that these are in fact the correct polynomials.

#figure(
  image("output/2k_lagrange_polynomials.png", width: 80%),
  caption: [Lagrange Polynomials $h_j (x)$ for $N=6$. Note that $h_j (x_i) = δ_(i j)$ where the abscisses $x_i$ are demarked by the dashed vertical lines.]
) <fig:lagrange_polynomials_2>

We then consider the duality between the nodal and the modal expansions for some function f,
$
bold(f)=bold(cal(V))hat(bold(f))
$
where $bold(f)=(f(tilde(x)_0),f(tilde(x)_1),...,f(tilde(x)_99))^top$ are the nodal values, and $hat(bold(f))$ denotes the modal values, see @eq:nodal_values. We can now approximate any function using the nodal values and the calculated Lagrange polynomials as our basis.

=== Approximating a Function Using Nodal Expansion
We do this for the function $v(x)=sin(pi x)$, $x in [-1,1]$. We first calculate the modal values as in h) and the Lagrange polynomials as above. Then we can calculate the interpolant,
$
I_N v(x)=sum_(j=0)^(N-1) (bold(cal(V))hat(bold(v)))_j thin h_j (x).
$

To inspect the result, we calculate the interpolant for varying $N$ and plot both the interpolant and the errors, see @fig:approximation_errors. Looking at the plots of the approximations, we see that we quite quickly get satisfying results (in the eyeball-norm). This is confirmed by the error plot, where we see that we reach errors of the size of the machine precision using as little as 21 grid points. We see a convergence that looks like the spectral convergence, we would expect for the DPT method, which makes sense, as the two expansions are equal.

#figure(
  image("output/2k_errors.png", width: 110%),
  caption: [Plot of the approximation of $v(x)$ for $x in [-1,1]$ for a selection of values of $N$ together with a plot of the truncation error $||v(x) -I_N v(x)||_2$ in a semi-log plot.]
) <fig:approximation_errors>

=== Moving Outside Interval
If we consider the expansion outside of the interval of our Jacobi polynomials, e.g. $x in [-1.5,1.5]$, the picture changes. We define the abscissas as before and do everything in the exact same way, except that we now do it for 150 uniform grid points $tilde(x)_i in [-1.5,1.5]$. We change the number of grid points to ensure that the distance between each pair is the same as before. 

We inspect the results like before, see @fig:approximation_errors_extended. Considering the plots of the interpolant, we see that errors occur outside $[-1,1]$ and that the interpolant to grow when moving further away. Looking at the plot of the error, we see that the error indeed grows rapidly for growing $N$ after $N=21$. 

By definition, the Jacobi polynomials only solve the Sturm-Liouville problem for $x in [-1,1]$ and since the abscissas, coefficients and Lagrange polynomials are based on these, it makes sense that the interpolant fails to be exact outside the interval.

#figure(
  image("output/2k_errors_extended.png", width: 110%),
  caption: [Plot of the approximation of $v(x)$ for $x in [-1.5,1.5]$ for a selection of values of $N$ together with a plot of the truncation error $||v(x) -I_N v(x)||_2$ in a semi-log plot.]
) <fig:approximation_errors_extended>

== Derivative of Jacobi Polynomials

In order to further leverage the power of the Vandermonde matrices, we employ a routine for finding the first derivatives of the Jacobi polynomials, which are given by @L2_slides[S.15]:
$
dv(,x) P_n^((α,β)) (x) = 1/2 (α + β + n + 1) P_(n-1)^((α+1, β+1)) (x)
$ <eq:jacobi_polynomials_deriv>

Where $P_n^((α,β))$ are the Jacobi Polynomials as given in @eq:jacobi_polynomials.
We implement a single routine which returns an $(n+1)×m$ matrix where $n$ is the order of the highest-order polynomial and $m$ is the number of nodes in $[-1, 1] ∋ x$ at which the derivatives are evaluated such. As such, each column represents a Jacobi polynomial evaluated at the nodes $x_i$.

This is implemented in Python as:

```python
def grad_jacobi_p(x: npt.NDArray, alpha: float, beta: float, n: int) -> npt.NDArray:
    """
    Computes the gradient of the first `n+1` Jacobi polynomials at nodes `x`.
    Reflects L2, slide 15.

    Arguments:
        x: Points at which to evaluate the gradients, shape (m,)
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$
        n: Highest order polynomial to compute (must be positive)
    Returns: Array of shape (m, n+1) where each column corresponds to the gradient of a Jacobi polynomial
    """
    grad_p = np.empty((len(x), n + 1))

    for i in range(n + 1):
        if i == 0:
            p_i = np.zeros_like(x)
        else:
            p_i = jacobi_p(x, alpha + 1, beta + 1, i - 1)[:, i - 1]

        coeff = 1 / 2 * (alpha + beta + i + 1)
        grad_p[:, i] = coeff * p_i

    return grad_p
```

Where the definition of the function `jacobi_p` has been relegated to @app:jacobi_p.

In order to gauge the correctness of our implementation, we plot the derivatives of the first 6 Legendre polynomials as shown in @fig:2k_grad_legendre.

#figure(
  image("output/2k_polynomial_grad_legendre.png"),
  caption: [
    Gradients of the first 6 Legendre polynomials evaluated over $x∈[-1, 1]$ using the `grad_jacobi_p` routine.
  ]
) <fig:2k_grad_legendre>

We may use this function to construct a generalized Vandermonde matrix $cal(V)_x$ enabling
the transformation between the nodal and modal coefficients, given by $f$ and  $hat(f)$ respectively, of a polynomial expansion of the first derivative of a function.

Using the same definition as in @sec:polynomial_j, we let $cal(V)_x$ be given by:
$
cal(V)_x ≔ mat(
  phi'_1(x_1), phi'_2(x_1), ..., phi'_N (x_1);
  phi'_1(x_2), phi'_2(x_2), ..., phi'_N (x_2);
  dots.v, dots.v, dots.down, dots.v;
  phi'_1(x_N), phi'_2(x_N), ..., phi'_N (x_N)
)
wider cal(V)_(x,i j) = eval(dv(φ_j, x))_(x=x_i)
$

We begin by considering the $N$-nodal expansion of the first derivative of $f(x)$ in the Lagrange basis, where $vv(f) = f(x_i)$ for $i ∈ 0, ..., N$:
$
dv(vv(f), x) ≈ dv(,x) (sum_(j=0)^N f_j h_j (x_i)) = sum_(j=0)^N f_j eval(dv(h_j,x) )_(x=x_i)
$

Where we introduce a _differentiation matrix_, $cal(D)$ characterised by its operation on the nodal coefficients $f_j$:
$
dv(vv(f), x)(x_i) = cal(D) vv(f)
$
$
cal(D)_(i j) ≔ eval(dv(h_j, x))_(x=x_i)
$ <eq:2k_diff_mat>

Problematically, the Lagrange polynomial construction is _global_ and as such the evaluation of $dv(h_j, x)$ at nodes $x_i$ requires information about all nodes, which may be computationally expensive. Instead, we would like to evaluate the derivative in the _modal_ basis, which is _local_ in the sense that a good modal representation will depend only on nearby 'neighbours', as seen through the recurrence relation of the Jacobi polynomials in @eq:jacobi_polynomials.

As such, we repeat the treatment above for the $N$-modal representation of the first derivative in a polynomial basis given by $φ_j (x)$:
$
dv(vv(f), x) ≈ dv(,x) (sum_(j=0)^N vv(hat(f))_j φ_j (x_i)) = sum_(j=0)^N vv(hat(f))_j eval(dv(φ_j,x) )_(x=x_i) = cal(V)_x vv(hat(f))
$ <eq:2k_modal>

As such, by recalling $vv(f) = cal(V) vv(hat(f))$ and using @eq:2k_diff_mat and @eq:2k_modal we find the relation:
$
dv(vv(f), x) cal(D) vv(f) = cal(D) cal(V) vv(hat(f)) = cal(V)_x vv(hat(f))
$
$
cal(D) = cal(V)_x cal(V)^(-1)
$

Where $cal(V)_x$ may be evaluated numerically using the function `grad_jacobi_p`.

We exercise our approach by solving the problem of evaluating the first derivative of function $v(x)$ given by:
$
v(x) = exp(sin(π x))
$ <eq:2k_v>

Where the domain of the function is given as $x ∈ [-1, 1]$. The analytical derivative of @eq:2k_v becomes:
$
dv(v, x) = π cos(π x) exp(sin(π x))
$

Using $N=$ nodes/modes we obtain the following following:
#figure(
  image("output/2k_grad_v.png"),
  caption: [
    First order derivative of $v(x) = sin(π x)$ evaluated on $x∈ [-1, 1]$ using the Vandermonde-derived differentiation matrix $cal(D)$ with $N=12$ nodes.
  ]
)

Lastly, we seek to evaluate the convergence of this method by computing the $L^2$-norm of the defect $dv(v, x)(x_i) - cal(D) vv(f)$ with $vv(f)$ being the nodal coefficients of $v(x_i)$.

This may be done by considering the definition of the $L^2$ norm:
$
norm(f)^2
&= ∫_Ω |f(x)|^2 dif x\
&≈ ∫_Ω vv(f)^TT vv(h)(x) (vv(f)^TT vv(h)(x))^TT dif x\
&= ∫_Ω vv(f)^TT vv(h)(x) vv(h)^TT (x) vv(f) dif x\
&= vv(f)^TT (∫_Ω (cal(V)^TT)^(-1) hat(vv(φ)) (x) hat(vv(φ))^TT (x) cal(V)^(-1) dif x) vv(f)\
&= vv(f)^TT (cal(V)^TT)^(-1) (∫_Ω hat(vv(φ)) (x) hat(vv(φ))^TT (x) dif x) cal(V)^(-1) vv(f)\
&= vv(f)^TT (cal(V)cal(V)^TT)^(-1) vv(f)\
$ <eq:2k_norm>

Where $Ω$ is the domain of the basis functions and the inner product of the basis functions evaluates to unity under the condition of orthonormality.

We may then define a _mass matrix_, $cal(M)$ as:
$
cal(M) ≔ (cal(V V^TT))^(-1)
$
Such that the $L^2$ norm may be evaluated as:
$
norm(f) ≈ sqrt(vv(f)^TT cal(M) vv(f))
$

We note that the function `jacobi_p` in @app:jacobi_p is not normalized, which is handled by implementing `jacobi_p_normalisation_const` and `jacobi_p_normalised` which may be found in @app:jacobi_p_normalized.

They correspond to the normalized Jacobi polynomials, $hat(P)_n^((α, β))$ @L2_slides[S.19]:
$
hat(P)_n^((α, β)) = 1/sqrt(γ_n^((α, β))) P_n^((α, β))
$

And associated normalisation constant $γ_n^((α, β))$ @L2_slides[S.11] given by:
$
γ_n^((α, β)) = 2^(α+β+1) (Γ(n + α + 1) Γ(n + β + 1))/(n! (2n + α + β + 1) Γ(n + α + β + 1))
$

Using these normalised polynomials to construct $cal(V)$ and subsequently $cal(M)$ allows us to easily evaluate the $L^2$-norm of the defects as a function of the number of modes, $N$, as shown in @fig:2k_convergence.

#figure(
  image("output/2k_convergence.png"),
  caption: [
    Convergence for approximation of first derivative of $v(x) = sin(π x)$ using differentiation matrix with Legendre polynomial basis in a semi-log plot.
  ]
) <fig:2k_convergence>

The expected convergence using a polynomial expansion will depend on how well-suited the polynomial series is to the modelled function, where particularly the _smoothness_ of the function may set a bound on the rate of convergence.

The Legendre polynomials leveraged here belong to the class of _ultraspherical polynomials_, $P_n^((α, α))$, while the function $v(x)$ belongs in the smoothness class $C^∞([-1, 1])$ by inspection. As such, borrowing the result from Hesthaven, Gottlieb & Gottlieb (2007) presented in @L2_slides[S.23-24] we find that we would expect to find _spectral convergence_, that is, the defect $d = dv(v, x)(x_1) - cal(D) vv(f)$, should follow:
$
norm(d) ≤ A e^(-B N)
$

Where $A, B ∈ ℝ_+$ are some constants dependant on the problem. This should form a straight line in a semi-log plot, as observed in @fig:2k_convergence, where we find a decent agreement with $B = 0.8$ as indicated by the dashed orange line.

We note that machine precision (≈ $10^(-13)$) is achieved with just $N=49$ modes.

We may conclude that the outlined method performs exceedingly well for sufficiently smooth functions, as evidenced by its performance on the sine function.


== Mass Matrix

In order to extend the domain $Ω$ of the norm $norm(f)$ from $[-1, 1]$ to an arbitrary interval $[a, b]$ where $a < b < ∞$,
we simply introduce an affine coordinate transformation into our construction of the mass matrix $cal(M)$ which scales and centers the domain:
$
macron(x)(x) = (b - a)/2 x + (a+b) /2
$ <eq:2l_transform>

Where $x$ is are the coordinates in the domain $Ω = [-1, 1]$ and $macron(x)$ the coordinates in the physical domain $macron(Ω) = [a, b]$.

We compute the differentials:
$
dif macron(x) = (b - a)/2 dif x
$

And again consider the definition of the $L^2$ inner product as done previously in @eq:2k_norm
and introduce a change of variables:
$
norm(f)^2
&= ∫_macron(Ω) |f(macron(x))|^2 dif macron(x)\
&= ∫_Ω |f(x)|^2 (b-a)/2 dif x \
&= (b-a)/2 ∫_Ω |f(x)|^2 dif x \
&= (b-a)/2 vv(f)^TT (cal(V)cal(V)^TT)^(-1) vv(f)\
$

As such, we simply need to transform the coordinates obtained from our root-finding routine into the physical domain using @eq:2l_transform and scale the mass matrix $cal(M)$ by a factor of $(b-a)/2$ in order to enable the $L^2$-norm to be found for any arbitrary domain $macron(Ω) = [a, b]$.

We employ this technique on the interval $macron(Ω) = [0, 2]$ and the following functions:
$
u_1(x) &≔ 1\
u_2(x) &≔ sin(x)
$

Whose squares we additionally integrate analytically to validate the correctness of our implementation:
$
∫_a^b |u_1(x)|^2 dif x = ∫_a^b 1 dif x = b - a
$
$
∫_a^b |u_2(x)|^2 dif x = ∫_a^b sin(x) dif x = [x/2 - sin(2x)/4]_a^b 
$

We again use the orthonormal Jacobi polynomials, $hat(P)_n^((α, β))$, to determine the Vandermonde matrix to ensure the correct scaling of the mass matrix.

Using $N$ nodes we find the following $L^2$-norms of $u_1 (x)$ and $u_2 (x)$ over $macron(Ω)=[0, 2]$:

#figure(
  table(
    columns: 5,
    table.header([*Function*], [*$N$*], [*Numerical*], [*Analytical*], [*Error*]),
    $u_1(x) = 1$, $5$, $1.4142135623730951$, $1.4142135623730951$, $<10^(-16)$,
    $u_2(x) = sin(x)$, $5$, $1.0905052416391159$, $1.0905047564439974$, $4.9 × 10^(-7)$,
    $u_2(x) = sin(x)$, $9$, $1.0905047564439974$, $1.0905047564439974$, $5.1 × 10^(-15)$,
  ),
  caption: [
    Comparison of numerical approximations of the $L^2$ norms of two functions, $u_1 (x)$, $u_2 (x)$
  ]
) <table:2l>

We note that in the definition of the mass matrix above and the derivation of its use in determining the $L^2$ norm, we introduce an approximation only in the second step of @eq:2k_norm:
$
norm(f)^2
&= ∫_Ω |f(x)|^2 dif x\
&≈ ∫_Ω vv(f)^TT vv(h)(x) (vv(f)^TT vv(h)(x))^TT dif x\
$

As such, our matrix-based integration using the mass matrix is exact under the assumption that the interpolating polynomials $vv(h)(x)$ (and equivalently, $hat(vv(φ))(x)$) are perfect representations of the modelled function. This will be the case if and only if the modelled function is itself a polynomial of order $q < p$ where $p = N-1$ is the highest order of the interpolating polynomials.

As such, we find that the $L^2$-norm of $u_1(x) = 1$ will be determined exactly using the numerical integration for interpolants with any number of nodes $N$, whereas $u_2(x) = sin(x)$ is only exact in the limit $N → ∞$, since the sine function is the sum of an infinite series of polynomials:
$
sin(x) = sum_(n=0)^∞ ((-1)^n x^(2n+1))/((2n+1)!)
$

Thankfully, as observed in @table:2l, we obtain machine precision at just $N=9$ nodes for $u_2(x)$.

#show: appendix
= Appendix

== Jacobi Polynomial Function <app:jacobi_p>

```python
def jacobi_p(x: npt.NDArray, alpha: float, beta: float, n: int) -> npt.NDArray:
    """
    Evaluates first `n+1` Jacobi polynomials at points `x` with parameters `alpha` and `beta`.
    Reflects L2, slide 12.

    Arguments:
        x: Points at which to evaluate the polynomials, shape (m,)
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$
        n: Highest order polynomial to compute (must be positive)

    Returns: Array of shape (m, n+1) where each column corresponds to a Jacobi polynomial
    """
    assert n >= 0, "n must be non-negative"

    P = np.empty((len(x), n + 1))

    P[:, 0] = 1.0
    if n == 0:
        return P

    P[:, 1] = 1 / 2 * (alpha - beta + (alpha + beta + 2) * x)
    if n == 1:
        return P

    for k in range(1, n):
        a_nm1_n = (
            2
            * (k + alpha)
            * (k + beta)
            / ((2 * k + alpha + beta + 1) * (2 * k + alpha + beta))
        )
        a_n_n = (alpha**2 - beta**2) / (
            (2 * k + alpha + beta + 2) * (2 * k + alpha + beta)
        )
        a_np1_n = (
            2
            * (k + 1)
            * (k + alpha + beta + 1)
            / ((2 * k + alpha + beta + 2) * (2 * k + alpha + beta + 1))
        )

        P[:, k + 1] = ((a_n_n + x) * P[:, k] - a_nm1_n * P[:, k - 1]) / a_np1_n

    return P
```

== Normalized Jacobi Polynomial Function <app:jacobi_p_normalized>
```python
def jacobi_p_normalisation_const(
    alpha: float, beta: float, n: int | npt.NDArray
) -> int | npt.NDArray:
    """
    Computes the normalisation constant for Jacobi polynomials.
    Reflects $γ_n^(α,β)$ from L2, slide 11.

    Arguments:
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$
        n: npt.ArrayLike

    Returns: normalisation constant(s) $γ$ for Jacobi polynomials $P_n$
    """
    return (
        2 ** (alpha + beta + 1)
        * (gamma(n + alpha + 1) * gamma(n + beta + 1))
        / (factorial(n) * (2 * n + alpha + beta + 1) * gamma(n + alpha + beta + 1))
    )


def jacobi_p_normalised(
    x: npt.NDArray, alpha: float, beta: float, n: int
) -> npt.NDArray:
    """
    Convenience function to get normalized Jacobi polynomials using `jacobi_p` and `jacobi_p_normalisation_const`.
    """
    P = jacobi_p(x, alpha, beta, n)
    norm_const = jacobi_p_normalisation_const(alpha, beta, np.arange(n + 1))
    return P / np.sqrt(norm_const)
```

#bibliography("references.bib") 