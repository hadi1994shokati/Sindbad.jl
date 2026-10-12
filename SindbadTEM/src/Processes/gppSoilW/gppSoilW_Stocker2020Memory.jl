export gppSoilW_Stocker2020Memory

#! format: off
@bounds @describe @units @timescale @with_kw struct gppSoilW_Stocker2020Memory{T1,T2,T3,T4} <: gppSoilW
    q::T1 = 1.0 | (0.01, 4.0) | "sensitivity of GPP to soil moisture " | "" | ""
    θstar::T2 = 0.6 | (0.1, 1.0) | "" | "" | ""
    τ_wmem_d::T3 = 365.0 | (30.0, 1500.0) | "recovery timescale of the demand-side drought memory (e-folding time of the stress memory)" | "day" | ""
    sens_wmem_d::T4 = 0.0 | (0.0, 1.0) | "sensitivity of the GPP soil-water factor to the demand-side drought memory; 0 = no memory (identical to gppSoilW_Stocker2020)" | "fraction" | ""
end
#! format: on

function define(params::gppSoilW_Stocker2020Memory, forcing, land, helpers)
    ## unpack parameters and land variables
    @unpack_gppSoilW_Stocker2020Memory params
    @unpack_nt z_zero ⇐ land.constants

    ## initial state: no remembered stress (equilibrated during spin-up)
    gpp_f_soilW = one(q)
    gpp_f_wmem_d = one(q)
    wmem_d = z_zero

    ## pack land variables
    @pack_nt wmem_d ⇒ land.states
    @pack_nt (gpp_f_soilW, gpp_f_wmem_d) ⇒ land.diagnostics
    return land
end

function compute(params::gppSoilW_Stocker2020Memory, forcing, land, helpers)
    ## unpack parameters
    @unpack_gppSoilW_Stocker2020Memory params

    ## unpack land variables
    @unpack_nt begin
        (∑w_fc, ∑w_wp) ⇐ land.properties
        soilW ⇐ land.pools
        wmem_d ⇐ land.states
        (z_zero, o_one, t_two) ⇐ land.constants
    end

    ## base factor: identical arithmetic to gppSoilW_Stocker2020
    SM = sum(soilW)
    max_AWC = at_least_zero(∑w_fc - ∑w_wp)
    actAWC = at_least_zero(SM - ∑w_wp)
    SM_nor = at_most_one(actAWC / max_AWC)
    tf_soilW = -q * (SM_nor - θstar)^t_two + o_one
    tmp_f_soilW = SM_nor <= θstar ? tf_soilW : one(tf_soilW)
    gpp_f_soilW_base = clamp_zero_one(tmp_f_soilW)

    ## drought memory: exponentially weighted mean of the daily soil-water stress (1 - normalised soil moisture)
    stress = o_one - SM_nor
    wmem_d = wmem_d + (stress - wmem_d) / τ_wmem_d

    ## legacy factor and combined soil-water factor (sens_wmem_d = 0 → factor 1 → base unchanged)
    gpp_f_wmem_d = o_one - sens_wmem_d * wmem_d
    gpp_f_soilW = gpp_f_soilW_base * gpp_f_wmem_d

    ## pack land variables
    @pack_nt wmem_d ⇒ land.states
    @pack_nt (gpp_f_soilW, gpp_f_wmem_d) ⇒ land.diagnostics
    return land
end

purpose(::Type{gppSoilW_Stocker2020Memory}) = "Soil moisture stress on GPP potential based on Stocker (2020), multiplied by a demand-side drought-memory factor driven by the remembered soil-water stress (test bed for learnable legacy effects)."

@doc """
$(getModelDocString(gppSoilW_Stocker2020Memory))

---

# Extended help

*Design*
 - Base factor: the exact arithmetic of `gppSoilW_Stocker2020` (`q`, `θstar`).
 - Memory state `wmem_d` (in `land.states`): exponentially weighted mean of the daily stress `1 - SM_nor`, where
   `SM_nor` is the normalised soil moisture of the whole column (WROASTED's actual `soilW` pools), with e-folding
   time `τ_wmem_d` (days; daily time step assumed). Initialised at 0 and equilibrated during spin-up.
 - Legacy factor `gpp_f_wmem_d = 1 - sens_wmem_d * wmem_d`, multiplied into `gpp_f_soilW`. With
   `sens_wmem_d = 0` the approach reproduces `gppSoilW_Stocker2020` bit for bit.
 - Companion supply-side hook: `transpirationSupply_wAWCMemory` (own state `wmem_s`, own parameters).

*Versions*
 - 1.0 on 12.10.2026 [DifferentiableHybridModel project]: W3 demand-side memory hook (H1)

*Created by*
 - DifferentiableHybridModel project (Hadi Shokati), with Claude
"""
gppSoilW_Stocker2020Memory
