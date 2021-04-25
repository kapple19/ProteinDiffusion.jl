H(x) = 2\(sign(x) + 1)

function findnearest(a,x)
	length(a) > 0 || return 0:-1
	r = searchsorted(a,x)
	length(r) > 0 && return r
	last(r) < 1 && return searchsorted(a,a[first(r)])
	first(r) > length(a) && return searchsorted(a,a[last(r)])
	x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
	x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
	return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end