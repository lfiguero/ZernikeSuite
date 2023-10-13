function normRatio()
	col = (α,j,n) -> [ZernikeSuite.ZernikePoly(α+j, i, n-i) for i=0:n];

	function ratio(f::ZernikeSuite.ZFun, j::Integer)
	    fshifted = lower(f)
	    for i = 2:j
		fshifted = lower(fshifted)
	    end
	    real(ZernikeSuite.wip(fshifted,fshifted)/ZernikeSuite.wip(f,f))
	end

	function randomguy(α::Real, j::Integer, n::Integer)
	    bf = col(α,j,n);
	    v = randn(n+1)+im*randn(n+1);
	    reduce(+, map(*, v, bf))
	end

	function randomratio(α::Real, j::Integer, n::Integer)
	    randomOP = randomguy(α,j,n);
	    ratio(randomOP, j)
	end

	maxratio = (α,j,n) -> maximum(map(g->ratio(g,j), col(α,j,n)));
	# So far I only know the formula for the j = 1 case
	myformula = (α,n) -> (n+1)/(α+1) + 1;
end
