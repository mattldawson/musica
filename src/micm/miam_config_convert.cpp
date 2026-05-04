// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Conversion from parsed MechanismConfiguration v1 aerosol types
// to MIAM ModelConfig types used by the MIAM builder.

#include <musica/micm/parse.hpp>

#include <mechanism_configuration/v1/aerosol_types.hpp>
#include <mechanism_configuration/v1/mechanism.hpp>

#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace musica
{
  namespace
  {
    namespace mc = mechanism_configuration::v1::types;

    template <typename T, typename = void>
    struct has_min_halflife : std::false_type
    {
    };

    template <typename T>
    struct has_min_halflife<T, std::void_t<decltype(std::declval<T>().min_halflife)>> : std::true_type
    {
    };

    /// @brief Find a species in the mechanism by name
    const mc::Species* FindSpecies(
        const std::string& name,
        const std::vector<mc::Species>& species)
    {
      for (const auto& s : species)
        if (s.name == name)
          return &s;
      return nullptr;
    }

    /// @brief Find a phase in the mechanism by name
    const mc::Phase* FindPhase(
        const std::string& name,
        const std::vector<mc::Phase>& phases)
    {
      for (const auto& p : phases)
        if (p.name == name)
          return &p;
      return nullptr;
    }

    /// @brief Find a PhaseSpecies in a phase by name
    const mc::PhaseSpecies* FindPhaseSpecies(
        const std::string& name,
        const mc::Phase& phase)
    {
      for (const auto& ps : phase.species)
        if (ps.name == name)
          return &ps;
      return nullptr;
    }

    miam_config::ArrheniusRateConstant ConvertArrhenius(const mc::AerosolArrheniusRate& rate)
    {
      return { rate.A, rate.C };
    }

    miam_config::EquilibriumConstant ConvertEquilibrium(const mc::AerosolEquilibriumConstant& eq)
    {
      return { eq.A, eq.C };
    }

    miam_config::HenrysLawConstant ConvertHLC(const mc::HenrysLawConstant& hlc)
    {
      return { hlc.hlc_ref, hlc.C };
    }
  }  // namespace

  std::optional<miam_config::ModelConfig> ConvertToMiamConfig(
      const mechanism_configuration::v1::types::Mechanism& mechanism)
  {
    if (mechanism.aerosol_representations.empty() && mechanism.aerosol_processes.empty())
      return std::nullopt;

    miam_config::ModelConfig config;
    config.name = mechanism.name;

    // ── Build species list ──────────────────────────────────────────
    for (const auto& spec : mechanism.species)
    {
      miam_config::SpeciesDef sd;
      sd.name = spec.name;
      sd.molecular_weight = spec.molecular_weight;
      config.species.push_back(sd);
    }

    // Also need to find density from phase-species and attach to config species
    // Build a map for quick lookup
    std::unordered_map<std::string, size_t> species_index;
    for (size_t i = 0; i < config.species.size(); ++i)
      species_index[config.species[i].name] = i;

    // ── Build phases ────────────────────────────────────────────────
    for (const auto& phase : mechanism.phases)
    {
      miam_config::PhaseDef pd;
      pd.name = phase.name;
      for (const auto& ps : phase.species)
      {
        pd.species_names.push_back(ps.name);

        // Propagate density from phase-species to config species
        if (ps.density.has_value())
        {
          auto it = species_index.find(ps.name);
          if (it != species_index.end())
            config.species[it->second].density = ps.density;
        }
      }

      if (phase.name == "gas")
        config.gas_phases.push_back(pd);
      else
        config.condensed_phases.push_back(pd);
    }

    // ── Convert representations ─────────────────────────────────────
    for (const auto& rep : mechanism.aerosol_representations)
    {
      std::visit(
          [&](const auto& r)
          {
            using T = std::decay_t<decltype(r)>;
            if constexpr (std::is_same_v<T, mc::UniformSection>)
            {
              miam_config::UniformSection us;
              us.name = r.name;
              us.phase_names = r.phases;
              us.min_radius = r.min_radius;
              us.max_radius = r.max_radius;
              config.representations.push_back(us);
            }
            else if constexpr (std::is_same_v<T, mc::SingleMomentMode>)
            {
              miam_config::SingleMomentMode sm;
              sm.name = r.name;
              sm.phase_names = r.phases;
              sm.geometric_mean_radius = r.geometric_mean_radius;
              sm.geometric_standard_deviation = r.geometric_standard_deviation;
              config.representations.push_back(sm);
            }
            else if constexpr (std::is_same_v<T, mc::TwoMomentMode>)
            {
              miam_config::TwoMomentMode tm;
              tm.name = r.name;
              tm.phase_names = r.phases;
              tm.geometric_standard_deviation = r.geometric_standard_deviation;
              config.representations.push_back(tm);
            }
          },
          rep);
    }

    // ── Convert processes and constraints ────────────────────────────
    for (const auto& proc : mechanism.aerosol_processes)
    {
      std::visit(
          [&](const auto& p)
          {
            using T = std::decay_t<decltype(p)>;

            if constexpr (std::is_same_v<T, mc::HenryLawEquilibrium>)
            {
              // HenryLawEquilibrium → HenryLawEquilibriumConstraint
              miam_config::HenryLawEquilibriumConstraint c;
              c.gas_species_name = p.gas_phase_species;
              c.condensed_species_name = p.condensed_phase_species;
              c.solvent_name = p.solvent;
              c.condensed_phase_name = p.condensed_phase;
              c.henrys_law_constant = ConvertHLC(p.henrys_law_constant);

              // Look up solvent MW and density from species/phase
              const mc::Species* solvent_spec = FindSpecies(p.solvent, mechanism.species);
              if (solvent_spec && solvent_spec->molecular_weight.has_value())
                c.mw_solvent = *solvent_spec->molecular_weight;
              else
                throw std::runtime_error(
                    "HENRY_LAW_EQUILIBRIUM: solvent '" + p.solvent + "' missing molecular weight");

              const mc::Phase* cond_phase = FindPhase(p.condensed_phase, mechanism.phases);
              if (cond_phase)
              {
                const mc::PhaseSpecies* solvent_ps = FindPhaseSpecies(p.solvent, *cond_phase);
                if (solvent_ps && solvent_ps->density.has_value())
                  c.rho_solvent = *solvent_ps->density;
                else
                  throw std::runtime_error(
                      "HENRY_LAW_EQUILIBRIUM: solvent '" + p.solvent +
                      "' missing density in phase '" + p.condensed_phase + "'");
              }

              config.constraints.push_back(c);
            }
            else if constexpr (std::is_same_v<T, mc::HenryLawPhaseTransfer>)
            {
              miam_config::HenryLawPhaseTransfer pt;
              pt.condensed_phase_name = p.condensed_phase;
              pt.gas_species_name = p.gas_phase_species;
              pt.condensed_species_name = p.condensed_phase_species;
              pt.solvent_name = p.solvent;
              pt.henrys_law_constant = ConvertHLC(p.henrys_law_constant);
              pt.diffusion_coefficient = p.diffusion_coefficient;
              pt.accommodation_coefficient = p.accommodation_coefficient;
              config.processes.push_back(pt);
            }
            else if constexpr (std::is_same_v<T, mc::DissolvedReaction>)
            {
              miam_config::DissolvedReaction dr;
              dr.phase_name = p.condensed_phase;
              dr.solvent_name = p.solvent;
              dr.rate_constant = ConvertArrhenius(p.rate_constant);
              if constexpr (has_min_halflife<T>::value)
                dr.min_halflife = p.min_halflife;
              for (const auto& r : p.reactants)
                dr.reactant_names.push_back(r.species_name);
              for (const auto& pr : p.products)
                dr.product_names.push_back(pr.species_name);
              config.processes.push_back(dr);
            }
            else if constexpr (std::is_same_v<T, mc::DissolvedReversibleReaction>)
            {
              miam_config::DissolvedReversibleReaction drr;
              drr.phase_name = p.condensed_phase;
              drr.solvent_name = p.solvent;
              for (const auto& r : p.reactants)
                drr.reactant_names.push_back(r.species_name);
              for (const auto& pr : p.products)
                drr.product_names.push_back(pr.species_name);
              if (p.forward_rate_constant)
                drr.forward_rate_constant = ConvertArrhenius(*p.forward_rate_constant);
              if (p.reverse_rate_constant)
                drr.reverse_rate_constant = ConvertArrhenius(*p.reverse_rate_constant);
              if (p.equilibrium_constant)
                drr.equilibrium_constant = ConvertEquilibrium(*p.equilibrium_constant);
              config.processes.push_back(drr);
            }
            else if constexpr (std::is_same_v<T, mc::DissolvedEquilibrium>)
            {
              miam_config::DissolvedEquilibriumConstraint dec;
              dec.phase_name = p.condensed_phase;
              dec.solvent_name = p.solvent;
              dec.algebraic_species_name = p.algebraic_species;
              dec.equilibrium_constant = ConvertEquilibrium(p.equilibrium_constant);
              for (const auto& r : p.reactants)
                dec.reactant_names.push_back(r.species_name);
              for (const auto& pr : p.products)
                dec.product_names.push_back(pr.species_name);
              config.constraints.push_back(dec);
            }
            else if constexpr (std::is_same_v<T, mc::LinearConstraint>)
            {
              miam_config::LinearConstraint lc;
              lc.algebraic_phase_name = p.algebraic_phase;
              lc.algebraic_species_name = p.algebraic_species;
              lc.constant = p.constant;
              lc.diagnose_from_state = p.diagnose_from_state;
              for (const auto& t : p.terms)
              {
                miam_config::LinearConstraintTerm term;
                term.phase_name = t.phase;
                term.species_name = t.species;
                term.coefficient = t.coefficient;
                lc.terms.push_back(term);
              }
              config.constraints.push_back(lc);
            }
          },
          proc);
    }

    return config;
  }

}  // namespace musica
