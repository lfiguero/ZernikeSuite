using PyPlot
include("rateAttainingSequences.jl")
switch_backend("agg") # We just want to export EPS files

# Run a small instance of the main program so Julia compiles it
conjecturedSharpnessTest(0.0, 2, 2:5);
list_of_α = [-0.99, 0.0, 9.9]
max_l = 3
number_of_sequence_members = 10
figcount = 0
colors = ["k", "b", "r", "g", "m"]
markers = ["o", "s", "^", "*", "D", "v", "h"]
cycledStyle(k::Int64) = colors[mod(k-1,length(colors))+1] * markers[mod(k-1,length(markers))+1] * "-"
cycledStyleAlt(k::Int64) = colors[mod(k-1,length(colors))+1] * "--"

for α in list_of_α
    for l in 1:max_l
	jarray = l + 2.^(1:number_of_sequence_members)
	ratios = conjecturedSharpnessTest(α, l, jarray)
	expectedRates = [-l, -1/2 + 2*(1:l) - l]
	N = 2jarray+2l-1
	# Plotting
	figure(figsize=(6,3))
	# Straight lines for rate comparison
	for k = 1:size(ratios,2)
	    # The C raises/lowers the straight lines showing the expected
	    # convergence/divergence rates so that they pass through their
	    # corresponding final datapoints
	    C = ratios[end,k] / N[end]^expectedRates[k]
	    loglog(N, C * N.^expectedRates[k], cycledStyleAlt(k), label="_nolegend_")
	end
	# Plots of 2j+2l-1 against the value of (residual norm-r / seminorm-l)
	# at the rate attaining sequence member of index j (logarithmic scale)
	for k = 1:size(ratios,2)
	    loglog(N, ratios[:,k], cycledStyle(k), label=LaTeXString("\$r = $(k-1)\$"))
	end
	legend(loc=0)
	title(LaTeXString("\$\\alpha = $α,\\ l = $l\$"))
	xlabel(LaTeXString("\$2j+2l-1\$"))
	ylabel("Seminorm ratio")
	tight_layout()
	# Save figure
	figcount = figcount + 1
	savefig("output/cst" * @sprintf("%03d", figcount) * ".eps")
	# Computation of empirical rates and export of part of a LaTeX tabular
	ER = zeros(size(ratios))
	ER[1,:] = nan(Float64)
	for k = 1:size(ratios, 2)
	    ER[2:end,k] = empiricalRate(N, ratios[:,k])
	end
	f = open("output/cst" * @sprintf("%03d", figcount) * ".tex", "w")
	table = """
	% The commands \\ratio and \\egr appearing below are not standard LaTeX
	% commands; it is up to the user to either replace or define them.
	"""
	table = table * "% \$\\alpha = $α\$, \$l = $l\$\n"
	table = table * "\$2j+2l-1\$ & "
	for k = 1:size(ratios, 2)
	    table = table * "\$\\ratio{\\alpha}{l}{$(k-1)}{j}\$ & "
	    table = table * "\$\\egr{\\alpha}{l}{$(k-1)}{j}\$" * (k<size(ratios, 2)?" & ":"\\\\\n")
	end
	for row = 1:size(ratios, 1)
	    table = table * string(N[row]) * " & "
	    for k = 1:size(ratios, 2)
		table = table * @sprintf("%5.2e", ratios[row,k]) * " & "
		table = table * @sprintf("%5.3f", ER[row,k]) * (k<size(ratios, 2)?" & ":"\\\\\n")
	    end
	end
	write(f, table)
	close(f)
    end
end
