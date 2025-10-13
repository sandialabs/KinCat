#include "KinCat.hpp"
#include "KinCat_CommandLineParser.hpp"
#include "KinCat_Main_Util.hpp"

#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"

int main(int argc, char *argv[]) {
  using real_type = KinCat::real_type;
  using ordinal_type = KinCat::ordinal_type;
  using site_type = KinCat::site_type;
  using device_type = KinCat::device_type;
  using host_device_type = KinCat::host_device_type;

  ///
  /// parse input
  ///

  /// later, we create a high level input including domain, dictionary and workflow
  std::string input_filename("input.json");
  ordinal_type verbose(0), verbose_parse(0), verbose_iterate(0);

  KinCat::CommandLineParser opts("KinCat Main");
  opts.set_option<std::string>("input", "Input dictionary file name", &input_filename);

  ///
  /// verbose level
  ///   0 - general function shows nothing,
  ///       showMe function shows brief summary of the object
  ///   1 - general function shows nothing (reserved for users' verbose output)
  ///       showMe function shows some details without lenghy output
  ///   2 - general function shows function begin and end
  ///       showMe function shows lengthy details
  ///   3 - general function shows internal workflow
  opts.set_option<ordinal_type>("verbose", "Verbosity level", &verbose);
  opts.set_option<ordinal_type>("verbose-parse", "Verbosity level in parsing step (not yet sorted)", &verbose_parse);
  opts.set_option<ordinal_type>("verbose-iterate", "Verbosity level in KMC iteration", &verbose_iterate);
  {
    const bool r_val = opts.parse(argc, argv);
    if (r_val)
      return 0;
  }

  Kokkos::initialize(argc, argv);
  do {
    ordinal_type r_val(0);

    ///
    /// show the used execution space
    ///
    {
#if defined(KOKKOS_ENABLE_CUDA)
      if (std::is_same<typename device_type::execution_space, Kokkos::Cuda>::value)
        std::cout << "-- Kokkos::Cuda is used\n";
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
      if (std::is_same<typename device_type::execution_space, Kokkos::OpenMP>::value)
        std::cout << "-- Kokkos::OpenMP is used\n";
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
      if (std::is_same<typename device_type::execution_space, Kokkos::Serial>::value)
        std::cout << "-- Kokkos::Serial is used\n";
#endif
    }

    ///
    /// input parser
    ///
    KinCat::Input in;
    in.parse(input_filename, verbose_parse);
    /// batch mode
    const ordinal_type n_samples = in._n_samples;

    ///
    /// create lattice
    ///
    KinCat::Lattice<device_type> lattice(n_samples, in._n_cells_x, in._n_cells_y, in._edge_vectors, in._basis_vectors,
                                         in._n_species);
    const ordinal_type search_distance = std::max(in._n_cells_interaction_x, in._n_cells_interaction_y);
    lattice.createNeighborPattern(in._symmetry_operations, in._site_coordinates, search_distance, verbose);
    if (in._site_init_type == "random") {
      if (in._site_random_fill_ratio_array.size() > 0) {
        /// variation of initial random configuration
        lattice.randomizeSites(in._site_random_fill_ratio_array, in._site_random_seed);
      } else {
        /// same initial configuration is used
        lattice.randomizeSites(in._site_random_fill_ratio, in._site_random_seed);
      }
    } else if (in._site_init_type == "file") {
      lattice.readinSites(in._sites, verbose_parse);
      /// restarting with multiple samples is tricky.
      /// the simulation should record 1) the last lattice configuration of each sample, 2) time, 3) different rates as
      /// each sample may use different rates overriden by initial ensemble input.
      /// Handling the issue 1) and 2) is small modification of existing code of KinCat_Main_Util.cpp
      /// To resolve issue 3), a user need to make sure that she/he uses the same rate configurations as is used or
      /// make sure that the simulation uses the same parameter for the runtime uncertainty.

    }
    KINCAT_CHECK_ERROR((in._n_cells_interaction_x * 2) + 1 > in._n_cells_domain_x,
                       "Error: the domain size in x needs to be over twice the interaction range in x");
    KINCAT_CHECK_ERROR((in._n_cells_interaction_y * 2) + 1 > in._n_cells_domain_y,
                       "Error: the domain size in y needs to be over twice the interaction range in y");
    lattice.setDomain(in._n_cells_x, in._n_cells_y);
    lattice.showMe(std::cout, "Lattice", verbose);

    ///
    /// create dictionary
    ///
    KinCat::value_type_2d_view<real_type, device_type> rates_all;

    if (in._process_override_array.size() > 0) {
      /// a different sample uses different rates
      rates_all = KinCat::value_type_2d_view<real_type, device_type>(KinCat::do_not_init_tag("rates all"), n_samples,
                                                                     in._rates.extent(1));
    } else {
      /// samples use the same rates
      rates_all = KinCat::value_type_2d_view<real_type, device_type>(KinCat::do_not_init_tag("rates all"), 1,
                                                                     in._rates.extent(1));
    }

    /// duplicate rates to rates_all for modification
    {
      /// Kokkkos MD RangePolicy is a convenient tool but it is still ambiguous how they map thread units
      /// to multi dimensional array. We prefer to use simple range policy.
      const auto rates = Kokkos::subview(in._rates, 0, Kokkos::ALL());
      Kokkos::parallel_for(
          Kokkos::RangePolicy<typename device_type::execution_space>(0, rates_all.span()),
          KOKKOS_LAMBDA(const ordinal_type ij) {
            const ordinal_type m = rates.extent(0), i = ij / m, j = ij % m;
            rates_all(i, j) = rates(j);
          });
      Kokkos::fence();
    }
    KinCat::ProcessDictionary<device_type> dictionary(in._variant_orderings, in._symmetry_orderings, in._configurations, in._configuration_sets, in._processints,
                                                      in._constraints, in._process_constraints, rates_all, in._process_symmetries);
    dictionary.showMe(std::cout, "Dictionary", verbose);
    dictionary.showMe(std::cout, "Instances details", in._processes, verbose);
    r_val = dictionary.validateData(in._processes, in._process_constraints);
    if (r_val > 0) {
      std::cerr << "Error: dictionary validateData returns with non-zero return value (" << r_val << ")\n";
      break;
    }

    if (in._instance_override_array.size() != 0 || in._process_override_array.size() != 0) {
      dictionary.overrideProcessRates(in._process_override_array, in._process_rates_override_array,
                                    in._instance_override_array, in._instance_rates_override_array, in._process_index, in._processes);
    }

    if (verbose_parse > 2) {
      std::cout << "Final Rates for all Samples \n";
      if (rates_all.extent(0) != in._n_samples) {
        std::cout << "Rates not varied between samples \n";
      }
      //std::cout << "rates_all dims (0,1) : " << rates_all.extent(0) << ", " << rates_all.extent(1) << '\n';
      for (ordinal_type j = 0, j_end = rates_all.extent(1); j < j_end; j++) {
        std::cout << "Instance " << j << ": ";
        for (ordinal_type i = 0, i_end = rates_all.extent(0); i < i_end; i++) {
          std::cout << rates_all(i, j) << ",   ";
        }
        std::cout << '\n';
      }
    }

    ///
    /// create counter
    ///
    KinCat::ProcessCounter<device_type> counter(1, in._processes.size(), dictionary._processints);
    counter.initialize(n_samples);

    ///
    /// main workflow
    ///
    KinCat::SolverBase<device_type> *solver =
        KinCat::SolverFactory<device_type>::createSolver(in._solver_type, verbose);

    solver->initialize(in._is_solver_random_number_variation, in._solver_random_seed, in._solver_random_pool_size,
                       in._n_kmc_steps_per_kernel_launch, in._n_cells_interaction_x, in._n_cells_interaction_y, lattice,
                       dictionary, counter, verbose);
    solver->showMe(std::cout, in._solver_type, verbose);

    /// time range
    KinCat::value_type_1d_view<real_type, device_type> t(KinCat::do_not_init_tag("t"), n_samples);
    if (in._t_sites.extent(0) == 0) { //restart file not read-in
      Kokkos::deep_copy(t, in._t_begin); // Use this if want to restart all simulations at same t_begin ([t_begin, t_end, t_step])
    }
    else {
      Kokkos::deep_copy(t,in._t_sites); // Use this to restart all samples at their previously reached times. 
    }
    const auto t_host = Kokkos::create_mirror_view_and_copy(host_device_type(), t);
    real_type t_global = in._t_begin;
    const real_type t_end = in._t_end;
    const real_type t_dt = in._t_dt;
    real_type t_step = t_dt > 0 ? std::floor(t_global / t_dt) * t_dt : 0;
    t_step += t_dt; // add now, since later adding now comes after advance, rather than before...

    const ordinal_type n_kmc_kernel_launches = in._n_kmc_kernel_launches;
    const ordinal_type n_kmc_steps_per_kernel_launch = in._n_kmc_steps_per_kernel_launch;
    ordinal_type epoch = 0;

    const real_type dump_interval = in._dump_interval;
    const real_type stats_interval = in._stats_interval;
    real_type dump_step = t_global + dump_interval; // initial data already outputted in dump/stats initialize functions
    real_type stats_step = t_global + stats_interval;

    /// snapshot
    bool dump_flag(true);
    if (in._dump_filename == "none" or in._dump_filename == "None") {
      dump_flag = false;
    }
    KinCat::Dump<device_type> dump(in._dump_filename, "dump-batch-site.json", lattice);
    if (dump_flag) {
      dump.showMe(std::cout, "Dump", verbose);
      dump.initialize(t_host, verbose);
    }

    bool stats_flag(true);
    if (in._stats_filename == "none" or in._stats_filename == "None") {
      stats_flag = false;
    }
    KinCat::Stats<device_type> stats(in._stats_filename, in._stats_list, lattice, in._processes, counter);
    if (stats_flag) {
      stats.showMe(std::cout, "Stats", verbose);
      stats.initialize(t_host, verbose);
    }
    if (!dump_flag && !stats_flag) {
      std::cout << "WARNING: No output files (dump, stats) specified! \n";
    }

    try {
      Kokkos::RangePolicy<typename device_type::execution_space> batch_range_policy(0, n_samples);
      for (ordinal_type iter = 0; iter < n_kmc_kernel_launches && t_global < t_end; ++iter) {
        const real_type t_prev(t_global);

        if (t_global >= t_step) { // Only advance t_step if previous kernel reached it, ensures that kernels always
                                  // stop at requested intervals
          t_step += t_dt;
          while ((t_global >= t_step) && t_dt != 0) { // simulation lept past next t_step, step again...
            t_step += t_dt;
          }
        }
        solver->advance(t, t_step, n_kmc_steps_per_kernel_launch, verbose);
        counter.syncToHost();
        Kokkos::Min<real_type> reducer_value(t_global);
        Kokkos::parallel_reduce(
            batch_range_policy,
            KOKKOS_LAMBDA(const ordinal_type sid, real_type &update) {
              const real_type t_cur = t(sid);
              update = (update < t_cur ? update : t_cur);
            },
            reducer_value);
        KINCAT_CHECK_ERROR(t_prev == t_global && t_global < t_end, "Error: time does not advance");

        Kokkos::deep_copy(t_host, t);

        /// output snapshots and advance t_dt for all domains
        if (dump_flag) {
          if (t_global >= dump_step) {
            dump.snapshot(t_host, verbose);
            dump_step += dump_interval;
          }
        }
        if (stats_flag) {
          if (t_global >= stats_step) {
            stats.snapshot(t_host, verbose);
            stats_step += stats_interval;
          }
        }

        if (t_global >= t_step) {
          if (verbose_iterate) {
            std::cout << " epoch = " << epoch << ", minimum time reached = " << t_global << "\n";
          }
          if (verbose_iterate > 1) {
            counter.showMe(std::cout, "Counter", in._processes, verbose_iterate);
          }
          epoch += 1;
        }
      }

      if (dump_flag) {
        dump.finalize(t_host, verbose);
      }
      if (stats_flag) {
        stats.finalize(t_host, verbose);
      }
      if (t_global < t_end) {
        std::cout << "WARNING: Not all samples reached the final time!!! \nCheck if kernel or step limited! ";
      }

    } catch (const std::exception &e) {
      std::cerr << "Error: exception is caught during kmc run\n" << e.what() << "\n";
      if (dump_flag) {
        dump.finalize(t_host, verbose);
      }
      if (stats_flag) {
        stats.finalize(t_host, verbose);
      }
    }
  } while (false);
  Kokkos::finalize();
}