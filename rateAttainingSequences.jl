include("ZernikeSuite.jl")
using ZernikeSuite

function pochhammer(x::Float64, k::Integer)
    out = 1.0
    for i = 0:k-1; out = out * (x+i); end
    out
end

relativeError(a, b) = 2*abs(a-b)/(abs(a)+abs(b))

empiricalRate(x, y) = log(y[2:end]./y[1:end-1]) ./ log(x[2:end]./x[1:end-1])

# The polynomial $t^{(\alpha,l)}_j$ of the paper
function t_alj(α::Real, l::Integer, j::Integer)
    coef = (-1)^l * exp(lgamma(α+j+1) - lgamma(j+1))^2 * exp(lgamma(α+2j+2) - lgamma(α+2j+l+2))
    t = coef * ZernikePoly(α, j, j)
    for k = l:-1:1
	coef = coef * k/(k-l-1) * (α+j+l-k+1)^2/(j+l-k+1)^2 * (α+2j+2l-2k+3)/(α+2j+2l-2k+1) * (α+2j+l-k+1)/(α+2j+2l-k+2)
	t = t + coef * ZernikePoly(α, j+l-k+1, j+l-k+1)
    end
    t
end

# ------------------------------------------------------------------------
# Computation of known ratios revealing sharpness of approximation results
# ------------------------------------------------------------------------
# Let us define $N^l_j$ as $2j+2l-1$ and $R^{\alpha,l}_j$ as
# $t^{\alpha,l}_j - \operatorname{Proj}^{\alpha}_{N^l_j}(t^{\alpha,l}_j)$.
# We compute the ratios between, on the one hand,
# $\| R^{\alpha,l}_j \|_{\mathrm{L}^2_{\rho^\alpha}(B^2)}^2$
# and
# $| R^{\alpha,l}_j |_{\mathrm{W}^{1,2}_{\rho^\alpha}(B^2)}^2$
# (standard seminorm) and, on the other hand,
# $| t^{\alpha,k}_j |_{\mathrm{W}^{l,2}_{\rho^\alpha}(B^2)}^2$
# (non-standard, albeit equivalent to the standard, seminorm involving the
# "$z$" and the "$z^*$" derivatives). We have explicit formulae for all three
# quantities which we can compare with what we obtain here by numerical means.
#
# The inputs are the quantities (in the notation used above) $\alpha$
# (character U+03B1 in the code below), $l$ and an array with the indices $j$
# to sample.
#
# The outputs are two arrays with the corresponding rate attainment revealing
# ratios of the computed norms and seminorms, and three arrays with the
# relative error between the computed norms and seminorms and what one obtains
# via the explicit formulae we have for them.
function knownSharpnessTest(α::Real, l::Integer, jarray::AbstractVector{Int64})
    nsamples = length(jarray)
    ratio_l2_norm_in_num = zeros(Float64, nsamples)
    ratio_h1_seminorm_in_num = zeros(Float64, nsamples)
    RE_res_sq_norm = zeros(Float64, nsamples)
    RE_res_sq_seminorm = zeros(Float64, nsamples)
    RE_t_sq_seminorm = zeros(Float64, nsamples)
    for (index, j) in enumerate(jarray)
	t = t_alj(α, l, j)
	res = t - proj(t, 2j+2l-1)
	# Norms and seminorms
	l2normressq = real(wip(res, res))
	w1seminormressq = real(wip(dx(res), dx(res)) + wip(dy(res), dy(res)))
	wlncseminormtsq = w_nc_sobolev_sq_sn(t, l)
	# Known explicit formulae for the norms and seminorms computed above
	formula1 = π*gamma(α+1)^2 / (2j+2l+α+1) * exp(2*lgamma(α+j+l+1) - 2*lgamma(j+l+1)) * exp(2*lgamma(α+2j+l+1) - 2*lgamma(α+2j+2l+1))
	formula2 = 4*π*gamma(α+1)^2 / (α+1) * exp(lgamma(α+j+l+1) - lgamma(j+l+1))^2 * ((α+2j+2l+1)^2*(j+l)*(α+j+l+1)) / pochhammer(α+2j+l+1, l+1)^2
	formula3 = 0.0
	for q = 0:l
	    formula3 = formula3 + exp(lgamma(α+j+l-q+1) + lgamma(α+j+q+1) - lgamma(j+l-q+1) - lgamma(j+q+1))
	end
	formula3 = π * gamma(α+1)^2 / (2j+l+α+1) * formula3
	display(index)
	ratio_l2_norm_in_num[index] = sqrt(l2normressq/wlncseminormtsq)
	ratio_h1_seminorm_in_num[index] = sqrt(w1seminormressq/wlncseminormtsq)
	RE_res_sq_norm[index] = relativeError(l2normressq, formula1)
	RE_res_sq_seminorm[index] = relativeError(w1seminormressq, formula2)
	RE_t_sq_seminorm[index] = relativeError(wlncseminormtsq , formula3)
    end
    (ratio_l2_norm_in_num, ratio_h1_seminorm_in_num, RE_res_sq_norm, RE_res_sq_seminorm, RE_t_sq_seminorm)
end

# ---------------------------------------------------------------------------
# Computation of unknown ratios suggesting sharpness of approximation results
# ---------------------------------------------------------------------------
# Let us define $N^l_j$ and $R^{\alpha,l}_j$ as above. For all integers $r$
# between $0$ and $l$ we compute the ratio between
# $| R^{\alpha,l}_j |_{\mathrm{W}^{r,2}_{\rho^\alpha}(B^2)}^2$
# and
# $| t^{\alpha,k}_j |_{\mathrm{W}^{l,2}_{\rho^\alpha}(B^2)}^2$
# in both cases using the standard seminorms.
#
# The inputs are the quantities (in the notation used above) $\alpha$
# (character U+03B1 in the code below), $l$ and an array with the indices $j$
# to sample.
#
# The output is a (number of indices $j$ requested)×(l+1) array with the
# corresponding ratios of the computed seminorms.
function conjecturedSharpnessTest(α::Real, l::Integer, jarray::AbstractVector{Int64})
    nsamples = length(jarray)
    ratios = zeros(Float64, nsamples, l+1)
    for (index, j) in enumerate(jarray)
	t = t_alj(α, l, j)
	res = t - proj(t, 2j+2l-1)
	# Seminorms
	wseminormsressq = all_w_sobolev_sq_sn(res, l)
	wlseminormtsq = w_sobolev_sq_sn(t, l)
	display(index)
	ratios[index,:] = sqrt(wseminormsressq/wlseminormtsq)
    end
    ratios
end
