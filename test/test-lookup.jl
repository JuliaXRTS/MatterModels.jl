using MatterModels
using Random

RNG = Xoshiro(137)

gaussian(x, y) = exp(-x^2 / 2 - y^2)
a_arr = collect(range(-rand(RNG), rand(RNG); length = 4)) # cols
b_arr = collect(range(-rand(RNG), rand(RNG); length = 5)) # rows
vals = gaussian.(a_arr', b_arr)
lookup_methods = [
    MatterModels.InterpolExtrapol(),
    MatterModels.InterpolEndValue(),
    MatterModels.NearestInput(),
    MatterModels.BelowInput(),
    MatterModels.AboveInput(),
]
GLT = MatterModels.GridLookupTable(a_arr, b_arr, vals)
@testset "$method" for method in lookup_methods
    for a in a_arr
        for b in b_arr
            test_val = MatterModels.lookup(method, a, b, GLT)
            groundtruth = gaussian(a, b)
            @test isapprox(groundtruth, test_val)
        end
    end
end
