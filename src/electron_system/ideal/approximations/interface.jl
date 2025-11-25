abstract type AbstractApproximation end
Base.broadcastable(approx::AbstractApproximation) = Ref(approx)

abstract type AbstractResponseApproximation <: AbstractApproximation end

# most naive implementation
struct NoApprox <: AbstractResponseApproximation end

# using the zero temperature limit (betabar to infty)
struct ZeroTemperatureApprox <: AbstractResponseApproximation end

# using the limit for non-degenerated electron gases (betabar << 1)
struct NonDegenerated <: AbstractResponseApproximation end

# using the limit for degenerated electron gases (betabar >> 1)
struct Degenerated <: AbstractResponseApproximation end
