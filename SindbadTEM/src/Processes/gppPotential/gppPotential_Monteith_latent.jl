export gppPotential_Monteith_latent
#=
DifferentiableHybridModel project (option 1 / "latent ecosystem state", approved by Hadi 2026-10-28):
Monteith potential GPP multiplied by a slowly varying, per-year latent factor exp(ls_y). The year of the current
time step comes from the forcing `f_year_idx` (1 .. 30 = index of the calendar year within the study period); a
non-integer index (e.g. the mean-seasonal-cycle spin-up forcing) or an index outside 1..30 gives factor 1.
With all ls_y = 0 the approach reproduces gppPotential_Monteith bit for bit. In the hybrid, the ls_y are learned as
free per-site values (embedding) with an L2 + smoothness penalty (see experiments/wroasted_fdk/src/LatentState.jl).
=#
@bounds @describe @units @timescale @with_kw struct gppPotential_Monteith_latent{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21,T22,T23,T24,T25,T26,T27,T28,T29,T30,T31} <: gppPotential
    εmax::T1 = 1.0 | (0.632, 3.2) | "Maximum Radiation Use Efficiency" | "gC/MJ" | ""
    ls_01::T2 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 1 (0 = no effect)" | "" | ""
    ls_02::T3 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 2 (0 = no effect)" | "" | ""
    ls_03::T4 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 3 (0 = no effect)" | "" | ""
    ls_04::T5 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 4 (0 = no effect)" | "" | ""
    ls_05::T6 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 5 (0 = no effect)" | "" | ""
    ls_06::T7 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 6 (0 = no effect)" | "" | ""
    ls_07::T8 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 7 (0 = no effect)" | "" | ""
    ls_08::T9 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 8 (0 = no effect)" | "" | ""
    ls_09::T10 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 9 (0 = no effect)" | "" | ""
    ls_10::T11 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 10 (0 = no effect)" | "" | ""
    ls_11::T12 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 11 (0 = no effect)" | "" | ""
    ls_12::T13 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 12 (0 = no effect)" | "" | ""
    ls_13::T14 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 13 (0 = no effect)" | "" | ""
    ls_14::T15 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 14 (0 = no effect)" | "" | ""
    ls_15::T16 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 15 (0 = no effect)" | "" | ""
    ls_16::T17 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 16 (0 = no effect)" | "" | ""
    ls_17::T18 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 17 (0 = no effect)" | "" | ""
    ls_18::T19 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 18 (0 = no effect)" | "" | ""
    ls_19::T20 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 19 (0 = no effect)" | "" | ""
    ls_20::T21 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 20 (0 = no effect)" | "" | ""
    ls_21::T22 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 21 (0 = no effect)" | "" | ""
    ls_22::T23 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 22 (0 = no effect)" | "" | ""
    ls_23::T24 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 23 (0 = no effect)" | "" | ""
    ls_24::T25 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 24 (0 = no effect)" | "" | ""
    ls_25::T26 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 25 (0 = no effect)" | "" | ""
    ls_26::T27 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 26 (0 = no effect)" | "" | ""
    ls_27::T28 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 27 (0 = no effect)" | "" | ""
    ls_28::T29 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 28 (0 = no effect)" | "" | ""
    ls_29::T30 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 29 (0 = no effect)" | "" | ""
    ls_30::T31 = 0.0 | (-0.7, 0.7) | "log of the latent capacity factor in study year 30 (0 = no effect)" | "" | ""
end

function define(params::gppPotential_Monteith_latent, forcing, land, helpers)
    @unpack_gppPotential_Monteith_latent params
    gpp_f_latent = one(εmax)
    @pack_nt gpp_f_latent ⇒ land.diagnostics
    return land
end

function compute(params::gppPotential_Monteith_latent, forcing, land, helpers)
    @unpack_gppPotential_Monteith_latent params
    @unpack_nt (f_PAR, f_year_idx) ⇐ forcing
    ls = (ls_01, ls_02, ls_03, ls_04, ls_05, ls_06, ls_07, ls_08, ls_09, ls_10, ls_11, ls_12, ls_13, ls_14, ls_15, ls_16, ls_17, ls_18, ls_19, ls_20, ls_21, ls_22, ls_23, ls_24, ls_25, ls_26, ls_27, ls_28, ls_29, ls_30)
    yi = isfinite(f_year_idx) ? round(Int, f_year_idx) : 0
    gpp_f_latent = (1 <= yi <= 30 && abs(f_year_idx - yi) < eltype(f_year_idx)(0.01)) ? exp(ls[yi]) : one(εmax)
    gpp_potential = εmax * gpp_f_latent * f_PAR
    @pack_nt (gpp_potential, gpp_f_latent) ⇒ land.diagnostics
    return land
end

purpose(::Type{gppPotential_Monteith_latent}) = "Potential GPP (Monteith radiation use efficiency) multiplied by a learnable per-year latent capacity factor exp(ls_y), selected by the forcing f_year_idx (latent ecosystem state test bed)."

@doc """
$(getModelDocString(gppPotential_Monteith_latent))

---

*Design*
 - `gpp_potential = εmax · exp(ls_y) · f_PAR`, where y = round(f_year_idx) if that is an integer in 1..30, else the factor is 1 (spin-up with mean-seasonal-cycle forcing, or days outside the study period).
 - `ls_y` default 0, bounds ±0.7 (factor 0.50 .. 2.0). All zero ⇒ identical to `gppPotential_Monteith`.
 - Diagnostic `gpp_f_latent` holds the factor of the current day.

*Versions*
 - 1.0 on 28.10.2026 [DifferentiableHybridModel project]: latent ecosystem state (option 1)

*Created by*
 - DifferentiableHybridModel project (Hadi Shokati), with Claude
"""
gppPotential_Monteith_latent
