import StatsBase as SB
import StatsBase: Histogram
import HypothesisTests
import Polynomials

## KS test of 2 distributions
function KS_test(Xdata, Tdata, bins)
    densX = DiscreteDensity(Xdata,bins)
    densT = DiscreteDensity(Tdata,bins)
    p_value = HypothesisTests.pvalue(HypothesisTests.ApproximateTwoSampleKSTest(densX, densT))
    return p_value
end

## discrete pdf
function DiscreteDensity(data, bins)
    N = length(data)
    histX = SB.fit(Histogram, data,  bins, closed=:left);
    density = histX.weights/N;
    return density;
end

## poly fit 
function poly_fit(X, y, reg_order=2, round_digits=4)
    fitted = Polynomials.fit(X, y, reg_order) |> p -> round.(Polynomials.coeffs(p), digits=4) |> Polynomials.Polynomial
    return fitted
end