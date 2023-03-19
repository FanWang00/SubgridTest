import StatsBase as SB
import StatsBase: Histogram
import HypothesisTests
import Polynomials


@doc raw"""  KS_test(Xdata, Tdata, bins)

Kolmogorov – Smirnov test for similarity of 2 empirical distributions 
"""
function KS_test(Xdata, Tdata, bins)

    densX = DiscreteDensity(Xdata,bins)
    densT = DiscreteDensity(Tdata,bins)
    p_value = HypothesisTests.pvalue(HypothesisTests.ApproximateTwoSampleKSTest(densX, densT))
    return p_value
end


@doc raw""" DiscreteDensity(data, bins)

Probability density estimation from data
"""
function DiscreteDensity(data, bins)
        """
        Empirical distributions  
        """
    N = length(data)
    histX = SB.fit(Histogram, data,  bins, closed=:left);
    density = histX.weights/N;
    return density;
end

@doc raw""" poly_fit(X, y, reg_order=2, round_digits=4)

Polynomial fit with X, y and regression order 
"""
function poly_fit(X, y, reg_order=2, round_digits=4)

    """
    polynomial fit with X, y and regression order 
    """

    fitted = Polynomials.fit(X, y, reg_order) |> p -> round.(Polynomials.coeffs(p), digits=4) |> Polynomials.Polynomial
    return fitted
end