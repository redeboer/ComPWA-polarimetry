
breakup(m, m1, m2) = sqrtKallenFact(m, m1, m2) / (2 * m)

function F²(l, p, p0, d)
    pR = p * d
    p0R = p0 * d
    l == 0 && return 1.0
    l == 1 && return (1 + p0R^2) / (1 + pR^2)
    l != 2 && error("l>2 cannot be")
    return (9 + 3p0R^2 + p0R^4) / (9 + 3pR^2 + pR^4)
end


abstract type Lineshape end

@with_kw struct BreitWignerMinL{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
BreitWignerMinL(pars::T; kw...) where {T} = BreitWignerMinL(; pars, kw...)
function (BW::BreitWignerMinL)(σ)
    dR, dΛc = 1.5, 5.0 # /GeV
    m, Γ₀ = BW.pars
    @unpack l, minL = BW
    @unpack m1, m2, mk, m0 = BW
    sqrtσ = sqrt(σ)
    p, p0 = breakup(sqrtσ, m1, m2), breakup(m, m1, m2)
    q, q0 = breakup(m0, sqrtσ, mk), breakup(m0, m, mk)
    Γ = Γ₀ * (p / p0)^(2l + 1) * m / sqrt(σ) * F²(l, p, p0, dR)
    1 / (m^2 - σ - 1im * m * Γ) * (p / p0)^l * (q / q0)^minL *
    sqrt(F²(l, p, p0, dR) * F²(minL, q, q0, dΛc))
end

# BuggBreitWignerMinL
@with_kw struct BuggBreitWignerMinL{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
BuggBreitWignerMinL(pars::T; kw...) where {
    T<:NamedTuple{X,Tuple{Float64,Float64}}} where {X} =
    BuggBreitWignerMinL(; pars=merge(pars, (γ=1.1,)), kw...)
#
function (BW::BuggBreitWignerMinL)(σ)
    σA = mK^2 - mπ^2 / 2
    m, Γ₀, γ = BW.pars
    @unpack m1, m2 = BW
    Γ = (σ - σA) / (m^2 - σA) * Γ₀ * exp(-γ * σ)
    1 / (m^2 - σ - 1im * m * Γ)
end

# Flatte1405
@with_kw struct Flatte1405{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
#
Flatte1405(pars::T; kw...) where {T} = Flatte1405(; pars, kw...)
function (BW::Flatte1405)(σ)
    m, Γ₀ = BW.pars
    @unpack m1, m2, m0, mk = BW
    sqrtσ = sqrt(σ)
    p, p0 = breakup(sqrtσ, m1, m2), breakup(m, mπ, mΣ)
    p′, p0′ = breakup(sqrtσ, mπ, mΣ), breakup(m, mπ, mΣ)
    Γ1 = Γ₀ * (p / p0) * m / sqrt(σ)
    Γ2 = Γ₀ * (p′ / p0′) * m / sqrt(σ)
    Γ = Γ1 + Γ2
    1 / (m^2 - σ - 1im * m * Γ)
end


function updatepars(BW::Lineshape, pars)
    fiels = fieldnames(typeof(BW))
    values = [getproperty(BW, f) for f in fiels]
    return typeof(BW)(; NamedTuple{fiels}(values)..., pars)
end



@recipe function f(BW::Lineshape)
    xv = range((BW.m1 + BW.m2)^2, (BW.m0 - BW.mk)^2, length=300)
    sqrtσ = sqrt(σ)
    intensity(σ) = abs2(BW(σ)) *
                   breakup(sqrtσ, BW.m1, BW.m2) *
                   breakup(BW.m0, sqrtσ, BW.mk) / sqrtσ
    yv = intensity.(xv)
    (xv, yv ./ sum(yv) .* length(yv))
end
