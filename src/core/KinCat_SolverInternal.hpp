#ifndef __KINCAT_SOLVER_INTERNAL_HPP__
#define __KINCAT_SOLVER_INTERNAL_HPP__

// clang-format off
#include "KinCat_Util.hpp"
#include "KinCat_Lattice.hpp"
#include "KinCat_ProcessDictionary.hpp"
// clang-format on

namespace KinCat {

namespace Impl {

template <typename T1, typename T2, typename CompareType>
KOKKOS_INLINE_FUNCTION static T1 *lower_bound(T1 *first, T1 *last, const T2 &val, CompareType compare) {
  T1 *it;
  ordinal_type step = 0, count = last - first;
  while (count > 0) {
    it = first;
    step = count / 2;
    it += step;
    if (compare(*it, val)) {
      first = ++it;
      count -= step + 1;
    } else {
      count = step;
    }
  }
  return first;
}

template <typename DT>
KOKKOS_INLINE_FUNCTION static void
findEvent(const ordinal_type sid, const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
          const ordinal_type k0, const ordinal_type k1, const real_type cell_specific_rate_to_search,
          ordinal_type &idx_event, ordinal_type &idx_pattern, ordinal_type &idx_variant, real_type &sum_rates_cell,
          const ordinal_type verbose) {
  sum_rates_cell = 0;

  using key_type = typename ProcessDictionary<DT>::key_type;
  key_type key;

  auto pattern = Kokkos::subview(lattice._pattern, 0, Kokkos::ALL(), Kokkos::ALL());
  const ordinal_type len_pattern = lattice._pattern.extent(0);
  const ordinal_type len_key = lattice._pattern.extent(1);

  ordinal_type idx_variant_local(-1);
  for (ordinal_type i = 0, iend = len_pattern; i < iend; ++i) {
    pattern.assign_data(&lattice._pattern(i, 0, 0));

    for (ordinal_type j = 0, jend = len_key; j < jend; ++j) {
      ordinal_type l0 = k0 + pattern(j, 0), l1 = k1 + pattern(j, 1), lb = pattern(j, 2);
      lattice.adjustPeriodicBoundary(l0, l1);

      ordinal_type lid;
      lattice.getSiteIndex(l0, l1, lb, lid);

      key[j] = lattice._sites(sid, lid);
    }

    ordinal_type idx_configuration(0), idx_first_event(0);
    const bool found_configuration = dictionary.searchConfiguration(key, idx_configuration, idx_variant_local);
    KINCAT_CHECK_ERROR(found_configuration != true, "Error: configuration key is not found in the dictionary");

    const bool found_event = dictionary.searchFirstEvent(idx_configuration, idx_first_event);
    if (found_event) {
      const ordinal_type idx_conf_a = dictionary._processints(idx_first_event, 0);
      for (ordinal_type j = idx_first_event, jend = dictionary._processints.extent(0); j < jend; ++j) {
        if (idx_conf_a != dictionary._processints(j, 0))
          break;

        const ordinal_type idx_process = dictionary._processints(j, 2);
        bool is_continue(true);
        {
          const ordinal_type n_allowed_symmetries = dictionary._process_symmetries(idx_process, 0);
          for (ordinal_type k = 0; k < n_allowed_symmetries; ++k) {
            const ordinal_type qid = dictionary._process_symmetries(idx_process, k + 1);
            const ordinal_type pid = lattice._symmetry2pattern(qid);
            if (i == pid) {
              is_continue = false;
              break;
            }
          }
        }
        if (is_continue)
          continue;

        const auto idx_variant_b = dictionary._constraints(j, idx_variant_local);
        if (idx_variant_b >= 0) {
          const real_type prev_sum_rates_cell = sum_rates_cell;
          const ordinal_type n_samples = dictionary._rates.extent(0);
          sum_rates_cell += dictionary._rates(sid % n_samples, j);
          if (cell_specific_rate_to_search >= prev_sum_rates_cell && cell_specific_rate_to_search < sum_rates_cell) {
            idx_event = j;
            idx_pattern = i;
            idx_variant = idx_variant_local;
          }
        }
      }
    }
  }
}

template <typename DT>
KOKKOS_INLINE_FUNCTION static void
computeCellRate(const ordinal_type sid, const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
                const ordinal_type k0, const ordinal_type k1, real_type &sum_rates_cell, const ordinal_type verbose) {
  const real_type cell_specific_rate_to_search(-1);
  ordinal_type idx_event, idx_pattern, idx_variant;
  findEvent(sid, lattice, dictionary, k0, k1, cell_specific_rate_to_search, idx_event, idx_pattern, idx_variant,
            sum_rates_cell, verbose);
}

template <typename DT>
KOKKOS_INLINE_FUNCTION static void
updateLatticeConfiguration(const ordinal_type sid, const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
                           const ordinal_type k0, const ordinal_type k1, const ordinal_type idx_event,
                           const ordinal_type idx_pattern, const ordinal_type idx_variant, const ordinal_type verbose) {
  const auto event = Kokkos::subview(dictionary._processints, idx_event, Kokkos::ALL());
  const auto conf_a = Kokkos::subview(dictionary._configurations, event(0), Kokkos::ALL());
  const auto conf_b = Kokkos::subview(dictionary._configurations, event(1), Kokkos::ALL());
  const auto pattern = Kokkos::subview(lattice._pattern, idx_pattern, Kokkos::ALL(), Kokkos::ALL());
  const auto idx_variant_b = dictionary._constraints(idx_event, idx_variant);

  for (ordinal_type j = 0, jend = pattern.extent(0); j < jend; ++j) {
    ordinal_type p0 = k0 + pattern(j, 0), p1 = k1 + pattern(j, 1), pb = pattern(j, 2);
    lattice.adjustPeriodicBoundary(p0, p1);

    ordinal_type pid;
    lattice.getSiteIndex(p0, p1, pb, pid);

    /// input double check
    {
      const ordinal_type k = dictionary._inverse_variant_orderings(idx_variant, j);
      KINCAT_CHECK_ERROR(lattice._sites(sid, pid) != conf_a(k), "Error: key (conf a) does not match to lattice");
    }

    /// update
    {
      const ordinal_type k = dictionary._inverse_variant_orderings(idx_variant_b, j);
      lattice._sites(sid, pid) = conf_b(k);
    }
  }
}

template <typename DT>
inline static void updateEventRatesLatticeDevice(const ordinal_type sid, const Lattice<DT> &lattice,
                                                 const ProcessDictionary<DT> &dictionary,
                                                 const value_type_1d_view<real_type, DT> &event_rates,
                                                 const ordinal_type verbose) {
  using range_policy_type = Kokkos::RangePolicy<typename DT::execution_space>;
  const range_policy_type policy(0, lattice.getNumberOfCells());

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ordinal_type cid) {
        ordinal_type k0, k1;
        lattice.getLatticeCellIndex(cid, k0, k1);
        lattice.adjustPeriodicBoundary(k0, k1);

        Impl::computeCellRate(sid, lattice, dictionary, k0, k1, event_rates(cid), verbose);
      });
}

template <typename DT>
inline static void updateEventRatesLatticeDevice(const ordinal_type sid, const Lattice<DT> &lattice,
                                                 const ProcessDictionary<DT> &dictionary,
                                                 const value_type_3d_view<real_type, DT> &event_rates,
                                                 const ordinal_type verbose) {
  using range_policy_type = Kokkos::RangePolicy<typename DT::execution_space>;
  const range_policy_type policy(0, lattice.getNumberOfCells());

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ordinal_type cid) {
        ordinal_type k0, k1;
        lattice.getLatticeCellIndex(cid, k0, k1);
        lattice.adjustPeriodicBoundary(k0, k1);

        ordinal_type d0, d1, lid;
        lattice.getDomainCellIndex(k0, k1, d0, d1, lid);

        Impl::computeCellRate(sid, lattice, dictionary, k0, k1, event_rates(d0, d1, lid), verbose);
      });
}

template <typename DT>
inline static void updateEventRatesLatticeDevice(const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
                                                 const value_type_2d_view<real_type, DT> &event_rates,
                                                 const ordinal_type verbose) {
  using range_policy_type = Kokkos::RangePolicy<typename DT::execution_space>;

  const ordinal_type n_samples = lattice.getNumberOfSamples();
  const ordinal_type n_cells = lattice.getNumberOfCells();
  const range_policy_type policy(0, n_samples * n_cells);

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ordinal_type iter) {
        const ordinal_type sid = iter / n_cells, cid = iter % n_cells;
        ordinal_type k0, k1;
        lattice.getLatticeCellIndex(cid, k0, k1);
        lattice.adjustPeriodicBoundary(k0, k1);

        Impl::computeCellRate(sid, lattice, dictionary, k0, k1, event_rates(sid, cid), verbose);
      });
}

template <typename DT>
KOKKOS_INLINE_FUNCTION static void
updateEventRates(const ordinal_type sid, const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
                 const ordinal_type arg_k0, const ordinal_type arg_k1, const ordinal_type n_cells_interaction_x,
                 const ordinal_type n_cells_interaction_y, const value_type_1d_view<real_type, DT> &event_rates,
                 const ordinal_type verbose) {
  /// range span
  const ordinal_type idx_begin_x = arg_k0 - n_cells_interaction_x + 1;
  const ordinal_type idx_end_x = arg_k0 + n_cells_interaction_x;
  const ordinal_type idx_begin_y = arg_k1 - n_cells_interaction_y + 1;
  const ordinal_type idx_end_y = arg_k1 + n_cells_interaction_y;

  /// for each cell
  for (ordinal_type k0_iter = idx_begin_x; k0_iter < idx_end_x; ++k0_iter) {
    for (ordinal_type k1_iter = idx_begin_y; k1_iter < idx_end_y; ++k1_iter) {
      ordinal_type k0 = k0_iter, k1 = k1_iter;
      lattice.adjustPeriodicBoundary(k0, k1);

      ordinal_type cid;
      lattice.getCellIndex(k0, k1, cid);

      Impl::computeCellRate(sid, lattice, dictionary, k0, k1, event_rates(cid), verbose);
    }
  }
}

template <typename DT>
KOKKOS_INLINE_FUNCTION static void
updateEventRates(const typename Kokkos::TeamPolicy<typename DT::execution_space>::member_type &member,
                 const ordinal_type sid, const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
                 const ordinal_type arg_k0, const ordinal_type arg_k1, const ordinal_type n_cells_interaction_x,
                 const ordinal_type n_cells_interaction_y, const value_type_3d_view<real_type, DT> &event_rates,
                 const ordinal_type verbose) {
  /// range span
  const ordinal_type idx_begin_x = arg_k0 - n_cells_interaction_x + 1;
  const ordinal_type idx_end_x = arg_k0 + n_cells_interaction_x;
  const ordinal_type idx_begin_y = arg_k1 - n_cells_interaction_y + 1;
  const ordinal_type idx_end_y = arg_k1 + n_cells_interaction_y;

  /// for each cell
  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, idx_begin_x, idx_end_x), [&](const ordinal_type &k0_iter) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, idx_begin_y, idx_end_y), [&](const ordinal_type &k1_iter) {
      ordinal_type k0 = k0_iter, k1 = k1_iter;
      lattice.adjustPeriodicBoundary(k0, k1);

      ordinal_type d0, d1, lid;
      lattice.getDomainCellIndex(k0, k1, d0, d1, lid);

      Impl::computeCellRate(sid, lattice, dictionary, k0, k1, event_rates(d0, d1, lid), verbose);
    });
  });
}

template <typename DT>
KOKKOS_INLINE_FUNCTION static void
updateEventRates(const ordinal_type sid, const Lattice<DT> &lattice, const ProcessDictionary<DT> &dictionary,
                 const ordinal_type arg_k0, const ordinal_type arg_k1, const ordinal_type n_cells_interaction_x,
                 const ordinal_type n_cells_interaction_y, const value_type_3d_view<real_type, DT> &event_rates,
                 const ordinal_type verbose) {
  /// range span
  const ordinal_type idx_begin_x = arg_k0 - n_cells_interaction_x + 1;
  const ordinal_type idx_end_x = arg_k0 + n_cells_interaction_x;
  const ordinal_type idx_begin_y = arg_k1 - n_cells_interaction_y + 1;
  const ordinal_type idx_end_y = arg_k1 + n_cells_interaction_y;

  /// for each cell
  for (ordinal_type k0_iter = idx_begin_x; k0_iter < idx_end_x; ++k0_iter) {
    for (ordinal_type k1_iter = idx_begin_y; k1_iter < idx_end_y; ++k1_iter) {
      ordinal_type k0 = k0_iter, k1 = k1_iter;
      lattice.adjustPeriodicBoundary(k0, k1);

      ordinal_type d0, d1, lid;
      lattice.getDomainCellIndex(k0, k1, d0, d1, lid);

      Impl::computeCellRate(sid, lattice, dictionary, k0, k1, event_rates(d0, d1, lid), verbose);
    }
  }
}

} // namespace Impl
} // namespace KinCat

#endif
