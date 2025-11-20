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

= Problem Formulation <sec:problem>

- Derive a SEM solver for the rotating hill problem. SEM requires a special type of discretization that is skew-symmetric for stability. This form is a result of a rewrite of the mathematical model using vector calculus, and then it is discretized.
- Solve the problem with a SEM in a circular domain (curved boundaries), i.e. a high-order scheme that has high-order and low-order in the same formulation.
- Perform convergence tests, i.e. h-convergence and p-convergence tests. 
- Use the solver to demonstrate break-even points for when high-order is more cost-efficient than low-order vs. high-order as a function of temporal integration time using an explicit 4th order runge-kutta method.

Questions for Allan:
- Governing equation (is there a document or paper we should use as reference?)
- Could he elaborate a bit on first point, hints on rewrite of mathematical model to get skew symmetry
- Hints/tools to do the mesh for a circular domain
  - Arises from the fact that rotation matrices are skew symmetric?
- What exactly is the difference between continuous and discontinuous Galerkin?
  - Intuition: how can we get $K N^p$ scaling if we impose continuous BC's between elements? Does this not couple the points of elements together?
- Should we use continuous or discontinuous method?

#bibliography("references.bib")
