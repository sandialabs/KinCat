#ifndef __KINCAT_TEST_RUNS_HPP__
#define __KINCAT_TEST_RUNS_HPP__

#include "KinCat_Main_Util.hpp"

TEST(Runs, CommonInput) {
  {
    const std::vector<std::string> input = {"test-files/input.json", "test-files/input-batch.json"};
    for (auto &filename : input) {
      std::cerr << ">> Testing filename : " << filename << "\n";
      KinCat::Input in;
      in.parse(filename);

      ordinal_type verbose(0), verbose_parse(0), verbose_iterate(0);
      ordinal_type n_samples(1);
      ordinal_type r_val(0);
      KinCat::Lattice<device_type> lattice(n_samples, in._n_cells_x, in._n_cells_y, in._edge_vectors, in._basis_vectors,
                                           in._n_species);
      const ordinal_type search_distance = std::max(in._n_cells_interaction_x, in._n_cells_interaction_y);
      lattice.createNeighborPattern(in._symmetry_operations, in._site_coordinates, search_distance, verbose);
      if (in._site_init_type == "random")
        lattice.randomizeSites(in._site_random_fill_ratio, in._site_random_seed);
      else if (in._site_init_type == "file")
        lattice.copySites(in._sites);

      lattice.setDomain(in._n_cells_domain_x, in._n_cells_domain_y);

      ///
      /// create dictionary
      ///
      KinCat::ProcessDictionary<device_type> dictionary(in._variant_orderings, in._configurations, in._processints,
                                                        in._constraints, in._rates, in._process_symmetries);
      dictionary.showMe(std::cout, "Dictionary", verbose);
      dictionary.showMe(std::cout, "Events details", in._processes, verbose);
      r_val = dictionary.validateData(in._processes, in._process_constraints);
      if (r_val > 0) {
        std::cerr << "Error: dictionary validateData returns with non-zero return value (" << r_val << ")\n";
        break;
      }

      ///
      /// create counter
      ///
      KinCat::ProcessCounter<device_type> counter(in._processes.size(), dictionary._processints);
      counter.initialize();

      ///
      /// main workflow
      ///
      KinCat::SolverBase<device_type> *solver = KinCat::SolverFactory<device_type>::createSolver(in._solver_type);

      solver->initialize(false, in._solver_random_seed, in._solver_random_pool_size,
                         in._n_kmc_steps_per_kernel_launch, in._n_cells_interaction_x, in._n_cells_interaction_y,
                         lattice, dictionary, counter, verbose);

      solver->showMe(std::cout, in._solver_type, verbose);

      /// time range
      ordinal_type n_domains, n_domains_x, n_domains_y;
      n_domains = lattice.getNumberOfDomains(n_domains_x, n_domains_y);
      KinCat::value_type_2d_view<real_type, device_type> t(KinCat::do_not_init_tag("t"), n_domains_x, n_domains_y);
      Kokkos::deep_copy(t, in._t_begin);

      KinCat::value_type_1d_view<real_type, host_device_type> t_global("t_global", n_samples);
      t_global(0) = in._t_begin;
      const real_type t_end = in._t_end;
      const real_type t_dt = in._t_dt;
      real_type t_step = t_dt > 0 ? std::floor(t_global(0) / t_dt) * t_dt : 0;

      const ordinal_type n_kmc_kernel_launches = in._n_kmc_kernel_launches;
      const ordinal_type n_kmc_iterations_per_kernel_launch = in._n_kmc_steps_per_kernel_launch;

      /// snapshot
      KinCat::Dump<device_type> dump(in._dump_filename, "dump-site.json", lattice);
      dump.showMe(std::cout, "Dump", verbose);

      dump.initialize(t_global, verbose);
      try {
        Kokkos::RangePolicy<typename device_type::execution_space> domain_range_policy(0, n_domains);
        for (ordinal_type iter = 0; iter < n_kmc_kernel_launches && t_global(0) < t_end; ++iter) {
          const real_type t_prev(t_global(0));
          /// advnace t_dt for all domains
          t_step += t_dt;
          solver->advance(t, t_step, n_kmc_iterations_per_kernel_launch, verbose);

          Kokkos::Min<real_type> reducer_value(t_global(0));
          Kokkos::parallel_reduce(
              domain_range_policy,
              KOKKOS_LAMBDA(const ordinal_type ij, real_type &update) {
                const ordinal_type i = ij / n_domains_y, j = ij % n_domains_y;
                const real_type t_cur = t(i, j);
                update = (update < t_cur ? update : t_cur);
              },
              reducer_value);
          KINCAT_CHECK_ERROR(t_prev == t_global(0) && t_global(0) < t_end, "Error: time does not advance");

          dump.snapshot(t_global, verbose);
          if (verbose_iterate) {
            std::cout << " iter = " << iter << ", n kmc kernel launches = " << n_kmc_iterations_per_kernel_launch
                      << ", " << n_kmc_kernel_launches << ", t = " << t_global(0) << ", t end = " << t_end << "\n";
          }
          if (verbose_iterate > 1)
            counter.showMe(std::cout, "Counter");
        }
        dump.finalize(t_global, verbose);
      } catch (const std::exception &e) {
        std::cerr << "Error: exception is caught during kmc run\n" << e.what() << "\n";
        dump.finalize(t_global, verbose);
      }
    }
    while (false)
      ;
  }
}

#endif
