"""
    metal_ff_absorption_departure(ν, T, Z, ni, ne, departure)

Computes the free-free linear absorption coefficient (in cm⁻¹) of a metal species for which a 
departure term is specified. (See [Peach1970](@ref) for more information.)

A free-free interaction is named as though the species interacting with the free electron had one 
more bound electron (in other words it's named as though the free-electron and ion were bound 
together). The `Z` argument and `ni` values should be specified for the species that actually
participates in the reaction. For completeness 2 examples provided below:
- Si I ff absorption: `ni` holds the number density of Si II, and `Z=1` (net charge of Si II)
- Si II ff absorption: `ni` holds the number density of Si III, and `Z=2` (net charge of Si III)

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `Z::Integer`: the net charge of the ion species (that participates in the interaction).
- `ni`: the number density of the ion species (that participates in the interaction) in cm⁻³.
- `ne`: the number density of free electrons.
- `departure`: a function or functor that gives the departure (i.e. relative difference from
  hydrogenic absorption) value at a given Temperature and σ, where 
  σ = (ν / Hz) / (Z² Rydberg_eV / hplanck_eV).
"""
function metal_ff_absorption_departure(ν::Real, T::Real, Z::Integer, ni::Real, ne::Real, departure)
    σ = ν/ Z^2 * (hplanck_eV / Rydberg_eV)
    hydrogenic_ff_absorption(ν, T, Z, ni, ne) * (1 + departure(T, σ))
end
#e.g.
#_He_I_ff(ν::Real, T::Real, ndens_HII::Real, nₑ::Real) =
#    metal_ff_absorption_departure(ν, T, 1, ndens_HII, nₑ, _departure_term_He_I_ff)


"""
    all_metal_ff_absorption(ν::Real, T::Real, number_densities::Dict, ne::Real,
                            exclude_species_with_departure_terms::Bool = false)

Computes the free-free linear absorption coefficient (in cm⁻¹) of metal species. This only includes
free-free absorption for species with a net positive charge of 1, 2, or 3. It also skips Hydrogenic
atoms.

# Arguments
- `ν`: frequency in Hz
- `T`: temperature in K
- `number_densities` is a `Dict` mapping each `Species` to its number density
- `ne`: the number density of free electrons.

This uses the hydrogenic approximation for all species except for those that have departure terms.

!!! note
    This mostly exists for demonstrative purposes.

    It might be useful to adopt this function for long-term use for efficiency reasons (even though
    it departs from the format of all other continuum absorption calculations).  

"""
function all_metal_ff_absorption(ν::Real, T::Real, number_densities::Dict, ne::Real,
                                 exclude_species_with_departure_terms::Bool = false)

    error("I think we may need the user to pass in the partition functions for every species")

    # recall, name of the free-free interaction differs from the name of the participating species
    departure_dict = Dict([(species"He_II", _He_I_ff),
                           (species"C_II",  _C_I_ff),
                           (species"C_III", _C_II_ff),
                           (species"Si_II", _Si_I_ff),
                           (species"Mg_II", _Mg_I_ff)])

    ndens_Z1 = 0.0
    ndens_Z2 = 0.0
    ndens_Z3 = 0.0

    α_out = 0.0 * ν

    for (k,ndens) in number_densities
        if k in [species"H_II", species"He_III", species"Li_IV"]
            # Honestly, it probably only makes sense to treat Hydrogen separately from this
            # function

            # skip cases where the species has not electrons, and is perfectly described by the
            # hydrogenic free-free absorption HI free-free, He II free-free, & Li III free-free
            continue
        elseif (exclude_species_with_departure_terms && (k in departure_dict))
            continue
        elseif k in departure_dict
            #add directly to α_out if there is a departure coefficient
            α_out += departure_dict[k](ν, T, ndens, ne)
        else
            #sum up contributions of hydrogenic ff coeffs, add them to α_out at the end
            if (k.charge == 1)     # example: species"O_II"
                ndens_Z1 += ndens
            elseif (k.charge == 2) # example: species"O_III"
                ndens_Z2 += ndens
            elseif (k.charge == 3) # example: species"O_IV"
                ndens_Z3 += ndens
            end
        end
    end

    #add contributions from species for which we use the uncorrected hydrogenic approximation
    α_out += hydrogenic_ff_absorption(T, ν, 1, ndens_Z1, ne)
    α_out += hydrogenic_ff_absorption(T, ν, 2, ndens_Z2, ne)
    α_out += hydrogenic_ff_absorption(T, ν, 3 #= matthew had this as 1=#, ndens_Z3, ne)

    α_out
end

