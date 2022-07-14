cd(joinpath(@__DIR__, ".."))
using Pkg
Pkg.activate(".")
Pkg.instantiate()
#
using YAML
#
using Parameters
#
using ThreeBodyDecay

using Lc2ppiKModelLHCb

using BenchmarkTools
using InteractiveUtils
using ThreadsX
using LinearAlgebra


import Plots.PlotMeasures.mm
using Plots
using LaTeXStrings

theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1, lab="", colorbar=false)

#                                  _|            _|
#  _|_|_|  _|_|      _|_|      _|_|_|    _|_|    _|
#  _|    _|    _|  _|    _|  _|    _|  _|_|_|_|  _|
#  _|    _|    _|  _|    _|  _|    _|  _|        _|
#  _|    _|    _|    _|_|      _|_|_|    _|_|_|  _|

isobarsinput = YAML.load_file(joinpath("..", "data", "particle-definitions.yaml"));

modelparameters =
    YAML.load_file(joinpath("..", "data", "model-definitions.yaml"));

const model0 = LHCbModel(
    modelparameters["Default amplitude model"];
    particledict=isobarsinput)

const λσs0 = randomPoint(tbs)
const σs0 = λσs0.σs
const two_λs = λσs0.two_λs


function intensity(m::LHCbModel, σs)
    @unpack chains, couplings = m
    #
    val = 0.0
    for two_λ0 in (-1, 1), two_λ1 in (-1, 1)
        two_λs = (two_λ1, 0, 0, two_λ0)
        av = amplitude.(chains, Ref(σs), Ref(two_λs))
        val += abs2(dot(av, couplings))
    end
    return val
end

intensity(model0, σs0)


# @code_warntype amplitude(model0.chains[1], σs0, two_λs)

# intensity(model0, σs0)


function benchmarkonset(model, Nev)
    ms = model.chains[1].tbs.ms
    data = flatDalitzPlotSample(ms; Nev)
    @belapsed intensity.(Ref($(model)), $(data))
end

function benchmarkonset_threads(model, Nev)
    ms = model.chains[1].tbs.ms
    data = flatDalitzPlotSample(ms; Nev)
    @belapsed ThreadsX.map(d -> intensity($(model), d), $(data))
end

sampleNev = [10, 100, 1000, 10000, 30000]
tv = benchmarkonset.(Ref(model0), sampleNev)
tv′ = benchmarkonset_threads.(Ref(model0), sampleNev)

let
    plot(xlab="Number of events", ylab="Computation time",
        ylim=(:auto, :auto), xscale=:log10, yscale=:log10)
    plot!(sampleNev, tv, lab="1 core", m=(:o, 3))
    plot!(sampleNev, tv′, lab="8 cores", m=(:o, 3))
end

# benchmarkonset_threads(model0, 100)


# benchmarkonset(model0, 100)
# benchmarkonset(model0, 1000)
# benchmarkonset(model0, 10000)
