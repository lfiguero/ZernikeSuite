using PyPlot
include("rateAttainingSequences.jl")
switch_backend("agg") # We just want to export EPS files

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
	figure(figsize=(4,3))
	plot(jarray, re1, "ko-", jarray, re2, "rs-", jarray, re3, "b^-", linewidth=1)
        legend(("RE 1", "RE 2", "RE 3"), loc=0)
	title(LaTeXString("\$\\alpha = " * string(α) * ",\\ l = " * string(l) * "\$"))
	xlabel(LaTeXString("\$j\$"))
	ylabel("Relative error")
	axis("tight")
	tight_layout()
	# Save figure
	figcount = figcount + 1
	savefig("output/kst" * @sprintf("%03d", figcount) * ".eps")
    end
end


