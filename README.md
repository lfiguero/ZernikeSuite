# ZernikeSuite: Utilities for manipulating finite linear combinations of Zernike polynomials

The main purpose of this repository is making available the numerical tests of the manuscript **Orthogonal polynomial projection error measured in Sobolev norms in the unit disk** <http://arxiv.org/abs/1503.04485v2>.

## Installation

This package should work for Julia 0.7 and above.

At the Julia prompt press `]` to enter the Pkg REPL mode and then enter:

```julia-repl
(v1.0) pkg> add https://github.com/lfiguero/ZernikeSuite.git
```

## `src/ZernikeSuite.jl`

In `src/ZernikeSuite.jl` lies the description of the ZFun type which represents functions as linear combinations of Zernike polynomials. A *tiny* number of operations are implemented—esentially, those used in the other codes mentioned below such as addition, some changes of parameters, projection, differentiation, some norms and seminorms, etc.

## `src/sharpness/sharpness.jl`

In `src/sharpness/sharpness.jl` lie the definitions of sequences of polynomials which depend on the parameters `α`, `l` and `j` (Equation (3.36) of the manuscript). They are meant to test the sharpness of a bound (Theorem 3.9 of the manuscript) on the (sometimes negative) growth with respect to `j` of a ratio of seminorms which depends on the same parameters `α`, `l` and `j` but also on a fourth parameter `r`. Then,

* if `r` is 0 or 1 there asymptotic growth rate of the seminorm ratio is known (that is, with a proof),
* if `r` is one of 2, 3, …, we have a conjectured asymptotic growth rate of the seminorm ratio but at this stage no proof for it.

In that same file the utilities in `src/ZernikeSuite.jl` are used to numerically compute those seminorm ratios and so:

* If `r` is 0 or 1 and
    + if one trusts the abovementioned proof, one can see the results as a validation of the numerical techniques; or
    + if one trusts the numerics, one can see the results as suggesting the proof is sound.
* On the other hand, if `r` is one of 2, 3, …, the numerics suggest that the conjectured asymptotic growth rate is sharp because is attained for the corresponding sequence of polynomials.

## Functions `runsConjecturedSharpnessTest` and `runsKnownSharpnessTest`

The functions `runsConjecturedSharpnessTest` and `runsKnownSharpnessTest` take as their argument the name of an output directory and perform test runs of the numerics of both knowledge regimes.

To run them, write in the Julia prompt, after installing this package,

```julia-repl
julia> using ZernikeSuite
julia> ZernikeSuite.runsConjecturedSharpnessTest("relative-path-of-an-output-directory")
julia> ZernikeSuite.runsKnownSharpnessTest("relative-path-of-an-output-directory")
```

The output consists of  plots and parts of LaTeX tabulars.
