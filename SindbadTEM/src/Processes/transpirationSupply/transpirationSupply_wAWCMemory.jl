export transpirationSupply_wAWCMemory

#! format: off
@bounds @describe @units @timescale @with_kw struct transpirationSupply_wAWCMemory{T1,T2,T3} <: transpirationSupply
    k_transpiration::T1 = 0.99 | (0.002, 1.0) | "fraction of total maximum available water that can be transpired" | "" | ""
    τ_wmem_s::T2 = 365.0 | (30.0, 1500.0) | "recovery timescale of the supply-side drought memory (e-folding time of the stress memory)" | "day" | ""
    sens_wmem_s::T3 = 0.0 | (0.0, 1.0) | "sensitivity of transpiration supply to the supply-side drought memory; 0 = no memory (identical to transpirationSupply_wAWC)" | "fraction" | ""
end
#! format: on

function define(params::transpirationSupply_wAWCMemory, forcing, land, helpers)
    @unpack_transpirationSupply_wAWCMemory params
    @unpack_nt z_zero ⇐ land.constants
    wmem_s = z_zero
    f_wmem_s = one(k_transpiration)
    transpiration_supply = z_zero
    @pack_nt wmem_s ⇒ land.states
    @pack_nt (transpiration_supply, f_wmem_s) ⇒ land.diagnostics
    return land
end

function compute(params::transpirationSupply_wAWCMemory, forcing, land, helpers)
    ## unpack parameters
    @unpack_transpirationSupply_wAWCMemory params

    ## unpack land variables
    @unpack_nt begin
        PAW ⇐ land.states
        (∑w_fc, ∑w_wp) ⇐ land.properties
        soilW ⇐ land.pools
        wmem_s ⇐ land.states
        o_one ⇐ land.constants
    end

    ## base supply: identical arithmetic to transpirationSupply_wAWC
    transpiration_supply_base = sum(PAW) * k_transpiration

    ## drought memory from the actual soil-water state (same stress definition as the demand-side hook)
    max_AWC = at_least_zero(∑w_fc - ∑w_wp)
    actAWC = at_least_zero(sum(soilW) - ∑w_wp)
    SM_nor = at_most_one(actAWC / max_AWC)
    stress = o_one - SM_nor
    wmem_s = wmem_s + (stress - wmem_s) / τ_wmem_s

    ## legacy factor on supply (sens_wmem_s = 0 → factor 1 → base unchanged)
    f_wmem_s = o_one - sens_wmem_s * wmem_s
    transpiration_supply = transpiration_supply_base * f_wmem_s

    ## pack land variables
    @pack_nt wmem_s ⇒ land.states
    @pack_nt (transpiration_supply, f_wmem_s) ⇒ land.diagnostics
    return land
end

purpose(::Type{transpirationSupply_wAWCMemory}) = "Supply-limited transpiration as in transpirationSupply_wAWC, multiplied by a supply-side drought-memory factor driven by the remembered soil-water stress (test bed for learnable legacy effects)."

@doc """
$(getModelDocString(transpirationSupply_wAWCMemory))

---

# Extended help

*Design*
 - Base supply: the exact arithmetic of `transpirationSupply_wAWC` (`k_transpiration`).
 - Memory state `wmem_s` (in `land.states`): exponentially weighted mean of the daily stress `1 - SM_nor`
   (normalised column soil moisture from WROASTED's `soilW` pools) with e-folding time `τ_wmem_s` (days).
 - Legacy factor `f_wmem_s = 1 - sens_wmem_s * wmem_s` multiplied into the supply. With `sens_wmem_s = 0`
   the approach reproduces `transpirationSupply_wAWC` bit for bit.
 - Companion demand-side hook: `gppSoilW_Stocker2020Memory` (own state `wmem_d`, own parameters).

*Versions*
 - 1.0 on 12.10.2026 [DifferentiableHybridModel project]: W3 supply-side memory hook (H2)

*Created by*
 - DifferentiableHybridModel project (Hadi Shokati), with Claude
"""
transpirationSupply_wAWCMemory
