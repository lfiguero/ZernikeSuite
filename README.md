# Utilities for manipulating finite linear combinations of Zernike polynomials

The main purpose of this repository is making available the numerical tests of the manuscript **Orthogonal polynomial projection error measured in Sobolev norms in the unit disk** <http://arxiv.org/abs/1503.04485v2>.

## ZernikeSuite.jl

In `ZernikeSuite.jl` lies the description of the ZFun type which represents functions as linear combinations of Zernike polynomials. A *tiny* number of operations are implemented—esentially, those used in the other codes mentioned below such as addition, some changes of parameters, projection, differentiation, some norms and seminorms, etc.

## rateAttainingSequences.jl

In `rateAttainingSequences.jl` lie the definitions of sequences of polynomials which depend on the parameters `α`, `l` and `j` (Equation (3.36) of the manuscript). They are meant to test the sharpness of a bound (Theorem 3.9 of the manuscript) on the (sometimes negative) growth with respect to `j` of a ratio of seminorms which depends on the same parameters `α`, `l` and `j` but also on a fourth parameter `r`. Then,

* if `r` is 0 or 1 there asymptotic growth rate of the seminorm ratio is known (that is, with a proof),
* if `r` is one of 2, 3, …, we have a conjectured asymptotic growth rate of the seminorm ratio but at this stage no proof for it.

In that same code the utilities in `ZernikeSuite.jl` are used to numerically compute those seminorm ratios and so:

* If `r` is 0 or 1 and
    + if one trusts the abovementioned proof, one can see the results as a validation of the numerical techniques; or
    + if one trusts the numerics, one can see the results as suggesting the proof is sound.
* On the other hand, if `r` is one of 2, 3, …, the numerics suggest that the conjectured asymptotic growth rate is sharp because is attained for the corresponding sequence of polynomials.

## runs-knownSharpnessTest.jl and runs-conjecturedSharpnessTest.jl

In `runs-knownSharpnessTest.jl` and `runs-conjecturedSharpnessTest.jl` there are test runs of the numerics of both knowledge regimes; when these codes are run (for example, by introducing `include(runs-knownSharpnessTest.jl` in Julia's interactive prompt) plots and parts of LaTeX tabulars with results are created in the `output` folder.

## How to run

Clone this repository (or use GitHub's *Download ZIP* button) and in the folder where the contents were stored run an interactive session of Julia and type `include("runs-knownSharpnessTest.jl")` or `include("runs-conjecturedSharpnessTest.jl")` (just one of them per session because their code import mechanisms collide).
Figures and pasteable portions of LaTeX tables are dumped in the `output` subfolder.
