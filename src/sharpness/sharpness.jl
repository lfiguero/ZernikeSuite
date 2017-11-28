##########################################################################
# Numerical experiments relevant for <http://arxiv.org/abs/1503.04485v2> #
##########################################################################
import PyPlot

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

function runsConjecturedSharpnessTest(outputDirectory::String)
	PyPlot.pygui(false) # We just want to export EPS files
	mkpath(outputDirectory)

	# Run a small instance of the main program so Julia compiles it
	conjecturedSharpnessTest(0.0, 2, 2:5);
	list_of_α = [-0.99, 0.0, 9.9]
	max_l = 3
	number_of_sequence_members = 12
	figcount = 0
	colors = ["k", "b", "r", "g", "m"]
	markers = ["o", "s", "^", "*", "D", "v", "h"]
	cycledStyle(k::Int64) = colors[mod(k-1,length(colors))+1] * markers[mod(k-1,length(markers))+1] * "-"
	cycledStyleAlt(k::Int64) = colors[mod(k-1,length(colors))+1] * "--"

	for α in list_of_α
	    for l in 1:max_l
		jarray = l + 2.^(1:number_of_sequence_members)
		ratios = conjecturedSharpnessTest(α, l, jarray)
		expectedRates = [-l; -1/2 + 2*(1:l) - l]
		Nlj = 2jarray+2l-1
		# Plotting
		PyPlot.figure(figsize=(6,3))
		# Straight lines for rate comparison
		for k = 1:size(ratios,2)
		    # The C raises/lowers the straight lines showing the expected
		    # convergence/divergence rates so that they pass through their
		    # corresponding final datapoints
		    C = ratios[end,k] / Nlj[end]^expectedRates[k]
		    PyPlot.loglog(Nlj, C * Nlj.^expectedRates[k], cycledStyleAlt(k), label="_nolegend_")
		end
		# Plots of 2j+2l-1 against the value of (residual norm-r / seminorm-l)
		# at the rate attaining sequence member of index j (logarithmic scale)
		for k = 1:size(ratios,2)
		    PyPlot.loglog(Nlj, ratios[:,k], cycledStyle(k), label=PyPlot.LaTeXString("\$r = $(k-1)\$"))
		end
		PyPlot.legend(loc=0)
		PyPlot.title(PyPlot.LaTeXString("\$\\alpha = $α,\\ l = $l\$"))
		PyPlot.xlabel(Pyplot.LaTeXString("\$N^{(l)}_j\$"))
		PyPlot.ylabel("Seminorm ratio")
		PyPlot.tight_layout()
		# Save figure
		figcount = figcount + 1
		PyPlot.savefig(outputDirectory * "/cst" * @sprintf("%03d", figcount) * ".eps")
		# Computation of empirical rates and export of part of a LaTeX tabular
		ER = zeros(size(ratios))
		ER[1,:] = convert(Float64, NaN)
		for k = 1:size(ratios, 2)
		    ER[2:end,k] = empiricalRate(Nlj, ratios[:,k])
		end
		f = open(outputDirectory * "/cst" * @sprintf("%03d", figcount) * ".tex", "w")
		table = """
		% The commands \\ratio, \\egr and \\nan appearing below are not
		% standard LaTeX commands; it is up to the user to either replace or
		% define them.
		"""
		table = table * "% \$\\alpha = $α\$, \$l = $l\$\n"
		table = table * "\$N^{(l)}_j\$ & "
		for k = 1:size(ratios, 2)
		    table = table * "\\ratio{\\alpha}{l}{$(k-1)}{j} & "
		    table = table * "\\egr{\\alpha}{l}{$(k-1)}{j}" * (k<size(ratios, 2)?" & ":"\\\\\n")
		end
		for row = 1:size(ratios, 1)
		    table = table * string(Nlj[row]) * " & "
		    for k = 1:size(ratios, 2)
			table = table * @sprintf("%5.2e", ratios[row,k]) * " & "
			table = table * @sprintf("%5.3f", ER[row,k]) * (k<size(ratios, 2)?" & ":"\\\\\n")
		    end
		end
		table = replace(table, "NaN", "\\nan")
		write(f, table)
		close(f)
	    end
	end
end

function runsKnownSharpnessTest(outputDirectory::String)
	PyPlot.pygui(false) # We just want to export EPS files
	mkpath(outputDirectory)

	# Run a small instance of the main program so Julia compiles it
	knownSharpnessTest(0.0, 2, 2:5);

	list_of_α = [-0.99, 0.0, 9.9]
	max_l = 3
	number_of_sequence_members = 100
	figcount = 0

	for i in 1:length(list_of_α)
	    for l in 1:max_l
		α = list_of_α[i]
		jarray = l:l+number_of_sequence_members-1
		ratio1, ratio2, re1, re2, re3 = knownSharpnessTest(α, l, jarray)
		expectedRate1 = -l
		expectedRate2 = 3/2 - l
		# Plots of relative errors between the numerical ratios and using the
		# known formulae
		PyPlot.figure(figsize=(4,3))
		PyPlot.plot(jarray, re1, "ko-", jarray, re2, "rs-", jarray, re3, "b^-", linewidth=1)
		PyPlot.legend(("RE 1", "RE 2", "RE 3"), loc=0)
		PyPlot.title(PyPlot.LaTeXString("\$\\alpha = " * string(α) * ",\\ l = " * string(l) * "\$"))
		PyPlot.xlabel(PyPlot.LaTeXString("\$j\$"))
		PyPlot.ylabel("Relative error")
		PyPlot.axis("tight")
		PyPlot.tight_layout()
		# Save figure
		figcount = figcount + 1
		PyPlot.savefig(outputDirectory * "/kst" * @sprintf("%03d", figcount) * ".eps")
	    end
	end
end
