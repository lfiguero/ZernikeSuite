# Utilities for manipulating finite linear combinations of Zernike polynomials

The main purpose of this repository is making available the numerical tests of a manuscript.

## ZernikeSuite.jl

In `ZernikeSuite.jl` lies the description of the ZFun type which represents functions as linear combinations of Zernike polynomials. A *tiny* number of operations are implemented—esentially, those used in the other codes mentioned below such as addition, some changes of parameters, projection, differentiation, some norms and seminorms, etc.

## rateAttainingSequences.jl

In `rateAttainingSequences.jl` lie the definition of a sequence of polynomials which depend on the parameters `α`, `l` and `j` which are meant to test the sharpness of a bound on the (sometimes negative) growth with respect to `j` of a ratio of seminorms which depends on the same parameters `α`, `l` and `j` but also on a fourth parameter `r`. Then,

* if `r` is 0 or 1 there asymptotic growth rate of the seminorm ratio is known (that is, with a proof),
* if `r` is one of 2, 3, …, we have a conjectured asymptotic growth rate of the seminorm ratio but at this stage no proof for it.

In that same code the utilities in `ZernikeSuite.jl` are used to numerically compute those seminorm ratios and so:

* If `r` is 0 or 1 and
    + if one trusts the abovementioned proof, one can see the results as a validation of the numerical techniques; or
    + if one trusts the numerics, one can see the results as suggesting the proof is sound.
* On the other hand, if `r` is one of 2, 3, …, the numerics suggest that the conjectured asymptotic growth rate is sharp because is attained for the corresponding sequence of polynomials.

## runs-knownSharpnessTest.jl and runs-conjecturedSharpnessTest.jl

In `runs-knownSharpnessTest.jl` and `runs-conjecturedSharpnessTest.jl` there are test runs of the numerics of both knowledge regimes; when these codes are run (for example, by introducing `include(runs-knownSharpnessTest.jl` in Julia's interactive prompt) plots and parts of LaTeX tabulars with results are created in the `output` folder.

