#include "KinCat.hpp"
#include "KinCat_CommandLineParser.hpp"
#include "KinCat_Main_Util.hpp"

int main(int argc, char *argv[]) {
  auto all_start = std::chrono::steady_clock::now();

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

    ///
    /// create lattice
    ///
    const ordinal_type n_samples(1);
    KinCat::Lattice<device_type> lattice(n_samples, in._n_cells_x, in._n_cells_y, in._edge_vectors, in._basis_vectors,
                                         in._n_species);
    const ordinal_type search_distance = std::max(in._n_cells_interaction_x, in._n_cells_interaction_y);
    lattice.createNeighborPattern(in._symmetry_operations, in._site_coordinates, search_distance, verbose);
    if (in._site_init_type == "random") {
      lattice.randomizeSites(in._site_random_fill_ratio, in._site_random_seed);
    } else if (in._site_init_type == "file") {
      lattice.readinSites(in._sites, verbose_parse);
    }

    KINCAT_CHECK_ERROR((in._n_cells_interaction_x * 2) + 1 > in._n_cells_domain_x,
                       "Error: the domain size in x needs to be over twice the interaction range in x");
    KINCAT_CHECK_ERROR((in._n_cells_interaction_y * 2) + 1 > in._n_cells_domain_y,
                       "Error: the domain size in y needs to be over twice the interaction range in y");
    
    lattice.setDomain(in._n_cells_domain_x, in._n_cells_domain_y);
    lattice.showMe(std::cout, "Lattice", verbose);

    ///
    /// create dictionary
    ///
    KinCat::ProcessDictionary<device_type> dictionary(in._variant_orderings, in._symmetry_orderings, in._configurations, in._configuration_sets, in._processints,
                                                      in._constraints, in._process_constraints, in._rates, in._process_symmetries);
    dictionary.showMe(std::cout, "Dictionary", verbose);
    dictionary.showMe(std::cout, "Instance details", in._processes, verbose);

    r_val = dictionary.validateData(in._processes, in._process_constraints);
    if (r_val > 0) {
      std::cerr << "Error: dictionary validateData returns with non-zero return value (" << r_val << ")\n";
      break;
    }

    ///
    /// create counter
    ///
    ordinal_type n_domains, n_domains_x, n_domains_y;
    n_domains = lattice.getNumberOfDomains(n_domains_x, n_domains_y);
    KinCat::ProcessCounter<device_type> counter(n_domains, in._processes.size(), dictionary._processints);
    counter.initialize();

    

    ///
    /// main workflow
    ///
    KinCat::SolverBase<device_type> *solver = KinCat::SolverFactory<device_type>::createSolver(in._solver_type);
    solver->initialize(false, in._solver_random_seed, in._solver_random_pool_size, in._n_kmc_steps_per_kernel_launch,
                       in._n_cells_interaction_x, in._n_cells_interaction_y, lattice, dictionary, counter, verbose);

    solver->showMe(std::cout, in._solver_type, verbose);

    /// time range
    KinCat::value_type_2d_view<real_type, device_type> t(KinCat::do_not_init_tag("t"), n_domains_x, n_domains_y);
    Kokkos::deep_copy(t, in._t_begin);

    KinCat::value_type_1d_view<real_type, host_device_type> t_global(
        "t_global", n_samples); // TODO: Should this be with domains, rather than samples? This is single sample version... Where is it used?

    t_global(0) = in._t_begin;
    const real_type t_end = in._t_end;
    const real_type t_dt = in._t_dt;
    real_type t_step = t_dt > 0 ? std::floor(t_global(0) / t_dt) * t_dt : 0; // limiting time for the next kernel

    const ordinal_type n_kmc_kernel_launches = in._n_kmc_kernel_launches;
    const ordinal_type n_kmc_steps_per_kernel_launch = in._n_kmc_steps_per_kernel_launch;
    const ordinal_type n_sites = lattice.getNumberOfSites();
    ordinal_type epoch = 0;

    const real_type dump_interval = in._dump_interval;
    const real_type stats_interval = in._stats_interval;
    real_type dump_step =
        t_global(0) + dump_interval; // initial data already outputted in dump/stats initialize functions
    real_type stats_step = t_global(0) + stats_interval;

    /// snapshot
    bool dump_flag(true);
    if (in._dump_filename == "none" or in._dump_filename == "None") {
      dump_flag = false;
    }
    #ifdef HAVE_HDF5
      bool useHDF5 = true;
    #else
    bool useHDF5 = false;
    #endif
    KinCat::Dump<device_type> dump(in._dump_filename, in._restart_save_filename, lattice, useHDF5);

    if (dump_flag) {
      dump.showMe(std::cout, "Dump", verbose);
      dump.initialize(t_global, verbose);
    }

    bool stats_flag(true);
    if (in._stats_filename == "none" or in._stats_filename == "None") {
      stats_flag = false;
    }
    KinCat::Stats<device_type> stats(in._stats_filename, in._stats_list, lattice, in._processes, counter, useHDF5);
    if (stats_flag) {
      stats.showMe(std::cout, "Stats", verbose);
      stats.initialize(t_global, verbose);
    }
    if (!dump_flag && !stats_flag) {
      std::cout << "WARNING: No output files (dump, stats) specified! \n";
    }

    ordinal_type old_proc_sum(0);
    real_type events_per_site(0.0);
    bool events_warning(false);

    try {
      Kokkos::RangePolicy<typename device_type::execution_space> domain_range_policy(0, n_domains);
      for (ordinal_type iter = 0; iter < n_kmc_kernel_launches && t_global(0) < t_end; ++iter) {
        const real_type t_prev(t_global(0));
        /// advance for all domains
        if (t_global(0) >= t_step) { // Only advance t_step if previous kernel reached it, ensures that kernels always
                                     // stop at requested intervals
          t_step += t_dt;
          while ((t_global(0) >= t_step) && t_dt != 0) { // simulation lept past next t_step, step again...
            t_step += t_dt;
          }
        }

        solver->advance(t, t_step, n_kmc_steps_per_kernel_launch, verbose);
        counter.syncToHost();
        //Check if sublattice algorithm has 'reasonable' steps
        if (in._solver_type == "sublattice") {
          ordinal_type proc_sum = counter.getTotalProcessCount(0);
          events_per_site = (real_type(proc_sum - old_proc_sum)/real_type(n_sites));
          old_proc_sum = proc_sum;
          if ((events_per_site >= 0.5)) {
            events_warning = true;
            std::cout << "Warning: averaged " << events_per_site << " events per site this kernel.\n";
          }
        }
        

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

        if (dump_flag) {
          if (t_global(0) >= dump_step) {
            dump.snapshot(t_global, verbose);
            dump_step += dump_interval;
          }
        }

        if (stats_flag) {
          if (t_global(0) >= stats_step) {
            stats.snapshot(t_global, verbose);
            stats_step += stats_interval;
          }
        }
	
        if (t_global(0) >= t_step) {
          if (verbose_iterate) {
            std::cout << " epoch = " << epoch << ", t = " << t_global(0) << "\n";
            counter.showMe(std::cout, "Counter", in._processes, verbose_iterate);
          }
          epoch += 1;
        }
      }

      if (dump_flag) {
        dump.finalize(t_global, verbose);
      }

      if (stats_flag) {
        stats.finalize(t_global, verbose);
      }

      if (t_global(0) < t_end) {
        std::cout << "WARNING: Simulation did not reach final time!!! \nCheck if kernel or step limited! \n";
      }
      if (events_warning) {
        std::cout << "WARNING: Simulation may have excessive errors due to sublattice solver settings. Carefully consider errors due to sublattice solver algorithm! Reduce kernel timestep to reduce errors.";
      }

    } catch (const std::exception &e) {
      std::cerr << "Error: exception is caught during kmc run\n" << e.what() << "\n";
      if (dump_flag) {
        dump.finalize(t_global, verbose);
      }

      if (stats_flag) {
        stats.finalize(t_global, verbose);
      }
    }
  } while (false);
  Kokkos::finalize();
}
