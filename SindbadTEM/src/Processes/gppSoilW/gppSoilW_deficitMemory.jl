export gppSoilW_deficitMemory

#! format: off
@bounds @describe @units @timescale @with_kw struct gppSoilW_deficitMemory{T1,T2,T3,T4,T5} <: gppSoilW
    τ_fast::T1 = 30.0 | (5.0, 120.0) | "timescale of the fast water-balance memory" | "day" | ""
    τ_slow::T2 = 365.0 | (60.0, 1500.0) | "timescale of the slow water-balance memory that carries the drought legacy" | "day" | ""
    wb_half::T3 = -1.0 | (-5.0, 2.0) | "water balance (rain - PET) at which the water limitation is 50 %" | "mm/day" | ""
    wb_slope::T4 = 2.0 | (0.5, 10.0) | "steepness of the water-limitation response" | "day/mm" | ""
    legacy_sens::T5 = 0.0 | (0.0, 1.0) | "sensitivity of GPP to the slow (legacy) water deficit; the default 0 switches the memory off" | "fraction" | ""
end
#! format: on

function define(params::gppSoilW_deficitMemory, forcing, land, helpers)
    ## unpack land variables
    @unpack_nt (z_zero, o_one) ⇐ land.constants

    ## initial memory states: neutral water balance, no limitation (equilibrated during spin-up)
    wb_fast = z_zero
    wb_slow = z_zero
    gpp_f_wb_fast = o_one
    gpp_f_wb_legacy = o_one
    gpp_f_soilW = o_one

    ## pack land variables
    @pack_nt (wb_fast, wb_slow) ⇒ land.states
    @pack_nt (gpp_f_soilW, gpp_f_wb_fast, gpp_f_wb_legacy) ⇒ land.diagnostics
    return land
end

function compute(params::gppSoilW_deficitMemory, forcing, land, helpers)
    ## unpack parameters
    @unpack_gppSoilW_deficitMemory params

    ## unpack forcing and land variables
    @unpack_nt f_rain ⇐ forcing
    @unpack_nt PET ⇐ land.fluxes
    @unpack_nt (wb_fast, wb_slow) ⇐ land.states
    @unpack_nt o_one ⇐ land.constants

    ## today's water balance (mm/day), negative = drying
    wb_now = f_rain - PET

    ## exponentially weighted memories of the water balance (daily time step; τ in days)
    wb_fast = wb_fast + (wb_now - wb_fast) / τ_fast
    wb_slow = wb_slow + (wb_now - wb_slow) / τ_slow

    ## current water limitation (0: full stress, 1: no stress) from the fast memory
    gpp_f_wb_fast = o_one / (o_one + exp(-wb_slope * (wb_fast - wb_half)))

    ## legacy limitation from the slow memory, scaled by the legacy sensitivity (in [1 - legacy_sens, 1])
    gpp_f_wb_legacy = o_one - legacy_sens * (o_one - o_one / (o_one + exp(-wb_slope * (wb_slow - wb_half))))

    ## combined effect passed on to gpp
    gpp_f_soilW = gpp_f_wb_fast * gpp_f_wb_legacy

    ## pack land variables
    @pack_nt (wb_fast, wb_slow) ⇒ land.states
    @pack_nt (gpp_f_soilW, gpp_f_wb_fast, gpp_f_wb_legacy) ⇒ land.diagnostics
    return land
end

purpose(::Type{gppSoilW_deficitMemory}) = "Water limitation of GPP with drought memory from fast and slow exponentially weighted water-balance (rain - PET) states; test bed for learning memory timescales (no soil water pools needed)."

@doc """
$(getModelDocString(gppSoilW_deficitMemory))

---

# Extended help

*Design*
 - Two states in `land.states` carry the past: `wb_fast` and `wb_slow` are exponentially weighted means of the
   daily water balance `f_rain - PET` with timescales `τ_fast` and `τ_slow` (days). They are initialised at 0
   (neutral) in `define` and equilibrate during spin-up.
 - `gpp_f_wb_fast = σ(wb_slope (wb_fast - wb_half))` is the current water limitation; `gpp_f_wb_legacy` is the
   same response applied to the slow state and scaled by `legacy_sens`, so that `legacy_sens = 0` switches the
   memory off. `gpp_f_soilW = gpp_f_wb_fast * gpp_f_wb_legacy` feeds the `gpp` approaches as usual.
 - Requires a `PET` approach earlier in the model structure (e.g. `PET_PriestleyTaylor1972`); intended for
   daily time steps (τ are in days and are not converted for other resolutions).
 - All operations are smooth, so the approach is differentiable with ForwardDiff through the time loop.

*Versions*
 - 1.0 on 09.10.2026 [DifferentiableHybridModel project]: first version, drought-memory test bed

*Created by*
 - DifferentiableHybridModel project (Hadi Shokati), with Claude

*Notes*
 - A test bed for the structure and learnability of memory, not a final process formulation.
"""
gppSoilW_deficitMemory
