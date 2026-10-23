# Write test/reference.json from the Julia version of the model, for test/model.test.mjs to check the JavaScript against.
# The Julia model lives in a separate, private repository; pass the path to its model.jl:
#     julia test/make_reference.jl /path/to/julia/model.jl
# The script needs no packages beyond the Julia standard library.

include(ARGS[1])

# The default set, and two others that move every parameter, so that a mistake in any one formula shows up.
const SETS = [
    PARAMS_FITTED,
    merge(PARAMS_FITTED, (fA = 0.05, alpha = 3.0, TT = 40.0, g1 = 0.3, kon2d = 0.5, kP = 0.4, nKP = 3.1,
                          lambdaX = 0.2, nC = 2.2, X_O1 = 0.03, O1max = 0.97)),
    PARAMS_NOMINAL,
]
const XS = 10.0 .^ (-8:0.5:2)
const KDS = [1.0, collect(values(PEPTIDE_KD))..., 1000.0]

json(x::Real) = isfinite(x) ? repr(Float64(x)) : "null"
json(x::AbstractVector) = "[" * join(json.(x), ", ") * "]"
json(nt::NamedTuple) = "{" * join(["\"$k\": $(json(v))" for (k, v) in pairs(nt)], ", ") * "}"

function reference(p)
    dp = discrimination_power(p)
    return (
        params = p,
        fraction = [cd69_fraction(x, KD, p) for KD in KDS for x in XS],
        log10_p15 = dp.log10_p15,
        slope = dp.alpha,
        intercept = dp.intercept,
        local_slopes = [local_p15_slope(KD, p) for KD in values(PEPTIDE_KD)],
    )
end

open(joinpath(@__DIR__, "reference.json"), "w") do io
    println(io, "{")
    println(io, "  \"xs\": ", json(XS), ",")
    println(io, "  \"KDs\": ", json(KDS), ",")
    println(io, "  \"sets\": [")
    println(io, join(["    " * json(reference(p)) for p in SETS], ",\n"))
    println(io, "  ]")
    println(io, "}")
end
