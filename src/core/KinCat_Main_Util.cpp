#include "KinCat_Main_Util.hpp"
#include "boost/json/src.hpp"
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace KinCat {

std::string replaceEnvironmentVariable(const std::string &in) {
  std::string r_val(in);
  std::regex env("\\$\\{([^}]+)\\}");
  std::smatch match;
  while (std::regex_search(r_val, match, env)) {
    const char *s = getenv(match[1].str().c_str());
    const std::string var(s == NULL ? "" : s);
    /// workaround for the compiler gcc 8 on silver
    const ordinal_type pos = match[0].first - r_val.begin();
    const ordinal_type len = match[0].second - match[0].first;
    r_val.replace(pos, len, var.c_str(), var.length());
  }
  return r_val;
}

ordinal_type parseInputfile(const std::string &filename, boost::json::value &jin, const ordinal_type verbose) {
  const FunctionScope func_scope("parseInputfile", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);
  try {
    std::ifstream file(filename);
    {
      std::string msg("Error: file is not open: " + filename);
      KINCAT_CHECK_ERROR(!file.is_open(), msg.c_str());
    }
    std::stringstream ss;
    ss << file.rdbuf();
    file.close();

    boost::json::error_code err;
    boost::json::parse_options option;
    option.allow_comments = true;
    option.allow_trailing_commas = true;
    jin = boost::json::parse(ss.str(), err, boost::json::storage_ptr(), option);
    {
      std::string msg("Error: boost json parsing error: " + err.message());
      KINCAT_CHECK_ERROR(err.value(), msg.c_str());
    }

    if (verbose > 3) {
      std::cout << "-- JSON after sanitizing\n";
      std::cout << jin << "\n";
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: exception is caught parsing JSON input file\n" << e.what() << "\n";
    r_val = -1;
  }

  return r_val;
}

ordinal_type validateDictionaryInput(const boost::json::value &jv, boost::json::object &root,
                                     const ordinal_type verbose) {
  const FunctionScope func_scope("validateDictionaryInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("dictionary");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }

    {
      const auto crystal_it = root.find("crystal");
      KINCAT_CHECK_ERROR(crystal_it == root.end(), "Error: [crystal] key does not exist in kincat input");
      {
        const auto crystal_obj = crystal_it->value().as_object();
        KINCAT_CHECK_ERROR(crystal_obj.find("edge vectors") == crystal_obj.end(),
                           "Error: [edge vectors] key does not exist in crystal object");
        KINCAT_CHECK_ERROR(crystal_obj.find("basis vectors") == crystal_obj.end(),
                           "Error: [basis vectors] key does not exist in crystal object");

        const auto symmetry_operation_it = crystal_obj.find("symmetry operations");
        KINCAT_CHECK_ERROR(symmetry_operation_it == crystal_obj.end(),
                           "Error: [symmetry operations] key does not exist in crystal object");
        {
          const auto symmetry_operation_obj = symmetry_operation_it->value().as_object();
          KINCAT_CHECK_ERROR(symmetry_operation_obj.find("shape") == symmetry_operation_obj.end(),
                             "Error: [shape] key does not exist in symmetry operation object");
          KINCAT_CHECK_ERROR(symmetry_operation_obj.find("data") == symmetry_operation_obj.end(),
                             "Error: [data] key does not exist in symmetry operation object");
        }
      }
    }

    {
      const auto configuration_it = root.find("configurations");
      KINCAT_CHECK_ERROR(configuration_it == root.end(), "Error: [configuration] key does not exist in kincat input");
      {
        const auto configuration_obj = configuration_it->value().as_object();
        KINCAT_CHECK_ERROR(configuration_obj.find("site coordinates") == configuration_obj.end(),
                           "Error: [site coordinates] key does not exist in configuration object");
        KINCAT_CHECK_ERROR(configuration_obj.find("variant orderings") == configuration_obj.end(),
                           "Error: [variant orderings] key does not exist in configuration object");
        KINCAT_CHECK_ERROR(configuration_obj.find("interaction range") == configuration_obj.end(),
                           "Error: [interaction range] key does not exist in configuration object");
        KINCAT_CHECK_ERROR(configuration_obj.find("shape") == configuration_obj.end(),
                           "Error: [shape] key does not exist in configuration object");
        KINCAT_CHECK_ERROR(configuration_obj.find("data") == configuration_obj.end(),
                           "Error: [data] key does not exist in configuration object");
        KINCAT_CHECK_ERROR(configuration_obj.find("configuration_sets") == configuration_obj.end(), 
                           "Error: [configuration_sets] key does not exist in configuration object")
      }
    }

    {
      const auto procdict_it = root.find("process dictionary");
      KINCAT_CHECK_ERROR(procdict_it == root.end(),
                         "Error: [process dictionary] key does not exist in dictionary input");
      {
        const auto procdict_obj = procdict_it->value().as_object();
        KINCAT_CHECK_ERROR(procdict_obj.find("processes") == procdict_obj.end(),
                           "Error: [processes] key does not exist in process dictionary object");
        KINCAT_CHECK_ERROR(procdict_obj.find("process symmetries") == procdict_obj.end(),
                           "Error: [process symmetries] key does not exist in process dictionary object");
        KINCAT_CHECK_ERROR(procdict_obj.find("process constraints") == procdict_obj.end(),
                           "Error: [process constraints] key does not exist in process dictionary object");
        KINCAT_CHECK_ERROR(procdict_obj.find("shape") == procdict_obj.end(),
                           "Error: [shape] key does not exist in process dictionary object");
        KINCAT_CHECK_ERROR(procdict_obj.find("data") == procdict_obj.end(),
                           "Error: [data] key does not exist in process dictionary object");
      }
    }

  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught validating input\n" << e.what() << "\n";
  }

  return r_val;
}

ordinal_type parseDictionaryInput(
    const boost::json::object &root, value_type_2d_view<real_type, device_type> &edge_vectors_view,
    value_type_2d_view<real_type, device_type> &basis_vectors_view,
    value_type_2d_view<real_type, device_type> &symmetry_operations_view,
    value_type_2d_view<real_type, device_type> &site_coordinates_view, ordinal_type &n_species,
    value_type_2d_view<ordinal_type, device_type> &variant_orderings_view, 
    value_type_2d_view<ordinal_type, device_type> &symmetry_orderings_view,
    ordinal_type &n_cells_interaction_x,
    ordinal_type &n_cells_interaction_y, value_type_2d_view<site_type, device_type> &configurations_view,
    value_type_1d_view<ordinal_type, device_type> &configuration_sets_view,
    value_type_2d_view<ordinal_type, device_type> &processints_view,
    value_type_2d_view<ordinal_type, device_type> &constraints_view,
    value_type_2d_view<ordinal_type, device_type> &process_symmetries_view, std::vector<std::string> &processes,
    std::map<std::string, ordinal_type> &process_index,
    std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints, const ordinal_type verbose) {
  const FunctionScope func_scope("parseDictionaryInput", __FILE__, __LINE__, verbose);

  std::vector<std::vector<real_type>> edge_vectors, basis_vectors;
  std::vector<ordinal_type> symmetry_operation_shape;
  std::vector<real_type> symmetry_operation_data;

  std::vector<std::vector<real_type>> site_coordinates;
  std::vector<std::vector<ordinal_type>> variant_orderings;
  std::vector<std::vector<ordinal_type>> symmetry_orderings;
  std::vector<ordinal_type> configuration_shape;
  std::vector<ordinal_type> configuration_data;
  std::vector<ordinal_type> configuration_sets;

  std::vector<std::vector<ordinal_type>> process_symmetries;
  processes.clear();
  process_constraints.clear();
  std::vector<ordinal_type> processints_shape;
  std::vector<ordinal_type> processints_data;

  try {
    ordinal_type r_val(0);
    const auto crystal_obj = root.find("crystal")->value().as_object();
    edge_vectors =
        boost::json::value_to<std::vector<std::vector<real_type>>>(crystal_obj.find("edge vectors")->value());
    basis_vectors =
        boost::json::value_to<std::vector<std::vector<real_type>>>(crystal_obj.find("basis vectors")->value());

    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "number of edge vectors : " << edge_vectors.size() << "\n";
      for (ordinal_type i = 0, iend = edge_vectors.size(); i < iend; ++i)
        std::cout << "  [" << edge_vectors[i][0] << "," << edge_vectors[i][1] << "]\n";
      std::cout << "\n";
      std::cout << "number of basis vectors : " << basis_vectors.size() << "\n";
      for (ordinal_type i = 0, iend = basis_vectors.size(); i < iend; ++i)
        std::cout << "  [" << basis_vectors[i][0] << "," << basis_vectors[i][1] << "]\n";
    }

    const auto symmetry_operation_obj = crystal_obj.find("symmetry operations")->value().as_object();
    r_val = KinCat::parseArrayObject(symmetry_operation_obj, symmetry_operation_shape, symmetry_operation_data);
    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "symmetry operations : (" << symmetry_operation_shape[0] << ", " << symmetry_operation_shape[1]
                << ")\n";
      std::cout << "  2x2 rotation matrix, 2x1 translation vector\n";
      for (ordinal_type i = 0, iend = symmetry_operation_shape[0]; i < iend; ++i) {
        const ordinal_type offs = i * symmetry_operation_shape[1];
        auto data = symmetry_operation_data.data() + offs;
        /// make it pretty output
        std::cout << "  rot : [ " << data[0] << ", " << data[1] << "; " << data[2] << ", " << data[3] << "], trans : ["
                  << data[4] << ", " << data[5] << "]\n";
      }
    }

    const auto configuration_obj = root.find("configurations")->value().as_object();
    site_coordinates =
        boost::json::value_to<std::vector<std::vector<real_type>>>(configuration_obj.find("site coordinates")->value());

    variant_orderings = boost::json::value_to<std::vector<std::vector<ordinal_type>>>(
        configuration_obj.find("variant orderings")->value());
    std::vector<ordinal_type> tmp =
        boost::json::value_to<std::vector<ordinal_type>>(configuration_obj.find("interaction range")->value());
    n_cells_interaction_x = tmp[0];
    n_cells_interaction_y = tmp[1];
    configuration_sets = boost::json::value_to<std::vector<ordinal_type>>(configuration_obj.find("configuration_sets")->value());
    ordinal_type n_sets = configuration_sets.size() - 1;
    std::vector<std::vector<ordinal_type>> symmetry_orderings_tmp;
    std::vector<ordinal_type> temp3;
    for (ordinal_type p = 0; p < variant_orderings.size(); p++) {
      symmetry_orderings.push_back(temp3);
    }
    const auto orderings_it = configuration_obj.find("symmetry orderings");
    if (orderings_it == configuration_obj.end()) {
      std::cout << "Cant find symmetry orderings\n";
    }
    const auto symmetry_orderings_obj = orderings_it->value().as_object();
    
    for (ordinal_type i = 0; i < n_sets; ++i) {
      symmetry_orderings_tmp = boost::json::value_to<std::vector<std::vector<ordinal_type>>>(
        symmetry_orderings_obj.find(std::to_string(i))->value());
      for (ordinal_type p = 0; p < symmetry_orderings_tmp.size(); p++) {
        for (ordinal_type j = 0; j < symmetry_orderings_tmp[0].size(); j++) {
          symmetry_orderings[p].push_back(symmetry_orderings_tmp[p][j]);
        }
      }
    }
    r_val = KinCat::parseArrayObject(configuration_obj, configuration_shape, configuration_data);
    
    {
      const auto it = std::max_element(configuration_data.begin(), configuration_data.end());
      n_species = *it + 1;
      KINCAT_CHECK_ERROR(n_species >= std::numeric_limits<site_type>::max(),
                         "Error: number of species is greater than the numeric limits of site type. site type can be "
                         "[char - 127], [short - 32767], [int - 2147483647]");
    }
    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "site coordinates : " << site_coordinates.size() << "\n";
      for (ordinal_type i = 0, iend = site_coordinates.size(); i < iend; ++i)
        std::cout << "  [" << site_coordinates[i][0] << ", " << site_coordinates[i][1] << "]\n";
      std::cout << "\n";
      std::cout << "interaction range : " << n_cells_interaction_x << ' ' << n_cells_interaction_y << "\n";
      std::cout << "variant orderings : " << variant_orderings.size() << "\n";
      for (ordinal_type i = 0, iend = variant_orderings.size(); i < iend; ++i) {
        std::cout << "  [" << variant_orderings[i][0];
        for (ordinal_type j = 1, jend = variant_orderings[i].size(); j < jend; ++j)
          std::cout << ", " << variant_orderings[i][j];
        std::cout << "]\n";
      }
      std::cout << "\n";
      std::cout << "symmetry orderings : " << symmetry_orderings.size() << "\n";
      for (ordinal_type i = 0, iend = symmetry_orderings.size(); i < iend; ++i) {
        std::cout << "  [" << symmetry_orderings[i][0];
        for (ordinal_type j = 1, jend = symmetry_orderings[i].size(); j < jend; ++j)
          std::cout << ", " << symmetry_orderings[i][j];
        std::cout << "]\n";
      }
      std::cout << "\n";
      std::cout << "configuration sets : [ " << configuration_sets[0];
      for (ordinal_type i = 1, iend = configuration_sets.size(); i < iend; i++) {
        std::cout << ", " << configuration_sets[i];
      }
      std::cout << " ]\n";
      std::cout << "configuration list : (" << configuration_shape[0] << ", " << configuration_shape[1] << ")\n";
      for (ordinal_type i = 0, iend = configuration_shape[0]; i < iend; ++i) {
        const ordinal_type offs = i * configuration_shape[1];
        std::cout << "  [" << i << ", [" << configuration_data[offs];
        for (ordinal_type j = 1, jend = configuration_shape[1]; j < jend; ++j)
          std::cout << ", " << configuration_data[offs + j];
        std::cout << "] ]\n";
      }
    }

    const auto procdict_obj = root.find("process dictionary")->value().as_object();
    processes = boost::json::value_to<std::vector<std::string>>(procdict_obj.find("processes")->value());
    process_index.clear();
    for (ordinal_type i = 0, iend = processes.size(); i < iend; ++i) {
      process_index[processes[i]] = i;
    }
    r_val = KinCat::parseArrayObject(procdict_obj, processints_shape, processints_data);
    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "processes : " << processes.size() << "\n";
      for (ordinal_type i = 0, iend = processes.size(); i < iend; ++i)
        std::cout << "  [" << i << ", \"" << processes[i] << "\"]\n";
      std::cout << "\n";
      std::cout << "instance list : (" << processints_shape[0] << ", " << processints_shape[1] << ")\n";
      std::cout << "  initial conf, final conf, process id\n";
      for (ordinal_type i = 0, iend = processints_shape[0]; i < iend; ++i) {
        const ordinal_type offs = i * processints_shape[1];
        std::cout << "  [" << i << ", [" << processints_data[offs + 0] << ", " << processints_data[offs + 1] << ", "
                  << processints_data[offs + 2] << "] ]\n";
      }
    }

    process_symmetries =
        boost::json::value_to<std::vector<std::vector<ordinal_type>>>(procdict_obj.find("process symmetries")->value());
    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "process symmetries : " << process_symmetries.size() << "\n";
      for (ordinal_type i = 0, iend = process_symmetries.size(); i < iend; ++i) {
        const auto &symmetries = process_symmetries[i];
        std::cout << "  [" << i << ", " << processes[i];
        { std::cout << ", [" << symmetries[0]; }
        for (ordinal_type j = 1, jend = symmetries.size(); j < jend; ++j) {
          std::cout << ", " << symmetries[j];
        }
        std::cout << "] ]\n";
      }
    }
    KINCAT_CHECK_ERROR(processes.size() != process_symmetries.size(),
                       "Error: processes size does not match to process symmetries");

    process_constraints = boost::json::value_to<std::vector<std::vector<std::vector<ordinal_type>>>>(
        procdict_obj.find("process constraints")->value());
    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "process constraints : " << process_constraints.size() << "\n";
      for (ordinal_type i = 0, iend = process_constraints.size(); i < iend; ++i) {
        const auto &constraints = process_constraints[i];
        std::cout << "  [" << i << ", " << processes[i];
        {
          const auto &constraint = constraints[0];
          std::cout << ", [ [" << constraint[0] << ", " << constraint[1] << ", " << constraint[2] << "]";
        }
        for (ordinal_type j = 1, jend = constraints.size(); j < jend; ++j) {
          const auto &constraint = constraints[j];
          std::cout << ", [" << constraint[0] << ", " << constraint[1] << ", " << constraint[2] << "]";
        }
        std::cout << " ] ]\n";
      }
    }

    KINCAT_CHECK_ERROR(processes.size() != process_constraints.size(),
                       "Error: processes size does not match to process constraints");

  } catch (const std::exception &e) {
    std::cerr << "Error: exception is caught parsing json input dictionary file\n" << e.what() << "\n";
  }

  /// create edge vectors view
  edge_vectors_view = value_type_2d_view<real_type, device_type>(do_not_init_tag("edge vectors"), edge_vectors.size(),
                                                                 edge_vectors[0].size());
  auto edge_vectors_host_view = Kokkos::create_mirror_view(host_device_type(), edge_vectors_view);
  copyStdVectorToKokkosView(edge_vectors, edge_vectors_host_view);
  Kokkos::deep_copy(edge_vectors_view, edge_vectors_host_view);

  /// basis vectors view
  basis_vectors_view = value_type_2d_view<real_type, device_type>(do_not_init_tag("basis vectors"),
                                                                  basis_vectors.size(), basis_vectors[0].size());
  auto basis_vectors_host_view = Kokkos::create_mirror_view(host_device_type(), basis_vectors_view);
  copyStdVectorToKokkosView(basis_vectors, basis_vectors_host_view);
  Kokkos::deep_copy(basis_vectors_view, basis_vectors_host_view);

  /// symmetry operations view
  symmetry_operations_view = value_type_2d_view<real_type, device_type>(
      do_not_init_tag("symmetry operations"), symmetry_operation_shape[0], symmetry_operation_shape[1]);
  auto symmetry_operations_host_view = Kokkos::create_mirror_view(host_device_type(), symmetry_operations_view);
  copyStdVectorToKokkosView(symmetry_operation_shape, symmetry_operation_data, symmetry_operations_host_view);
  Kokkos::deep_copy(symmetry_operations_view, symmetry_operations_host_view);

  /// site coordinates view
  site_coordinates_view = value_type_2d_view<real_type, device_type>(
      do_not_init_tag("site coordinates"), site_coordinates.size(), site_coordinates[0].size());
  auto site_coordinates_host_view = Kokkos::create_mirror_view(host_device_type(), site_coordinates_view);
  copyStdVectorToKokkosView(site_coordinates, site_coordinates_host_view);
  Kokkos::deep_copy(site_coordinates_view, site_coordinates_host_view);

  /// variant orderings view
  variant_orderings_view = value_type_2d_view<ordinal_type, device_type>(
      do_not_init_tag("variant orderings"), variant_orderings.size(), variant_orderings[0].size());
  auto variant_orderings_host_view = Kokkos::create_mirror_view(host_device_type(), variant_orderings_view);
  copyStdVectorToKokkosView(variant_orderings, variant_orderings_host_view);
  Kokkos::deep_copy(variant_orderings_view, variant_orderings_host_view);

  /// symmetry orderings view
  symmetry_orderings_view = value_type_2d_view<ordinal_type, device_type>(
    do_not_init_tag("symmetry orderings"), symmetry_orderings.size(), symmetry_orderings[0].size()); 
  auto symmetry_orderings_host_view = Kokkos::create_mirror_view(host_device_type(), symmetry_orderings_view);
  copyStdVectorToKokkosView(symmetry_orderings, symmetry_orderings_host_view);
  Kokkos::deep_copy(symmetry_orderings_view, symmetry_orderings_host_view);

  /// configuration and instance view
  configurations_view = value_type_2d_view<site_type, device_type>(do_not_init_tag("configurations"),
                                                                   configuration_shape[0], configuration_shape[1]);
  auto configurations_host_view = Kokkos::create_mirror_view(host_device_type(), configurations_view);
  copyStdVectorToKokkosView(configuration_shape, configuration_data, configurations_host_view);
  Kokkos::deep_copy(configurations_view, configurations_host_view);
  
  configuration_sets_view = value_type_1d_view<ordinal_type, device_type>(do_not_init_tag("configuration_sets"), configuration_sets.size());
  auto configuration_sets_host_view = Kokkos::create_mirror_view(host_device_type(), configuration_sets_view);
  copyStdVectorToKokkosView(configuration_sets, configuration_sets_host_view);
  Kokkos::deep_copy(configuration_sets_view, configuration_sets_host_view);

  processints_view = value_type_2d_view<ordinal_type, device_type>(do_not_init_tag("instances"), processints_shape[0],
                                                                   processints_shape[1]);
  auto processints_host_view = Kokkos::create_mirror_view(host_device_type(), processints_view);
  copyStdVectorToKokkosView(processints_shape, processints_data, processints_host_view);
  Kokkos::deep_copy(processints_view, processints_host_view);

  ordinal_type max_symmetries(0);
  for (const auto &s : process_symmetries) {
    const ordinal_type s_size = s.size();
    max_symmetries = max_symmetries > s_size ? max_symmetries : s_size;
  }
  process_symmetries_view = value_type_2d_view<ordinal_type, device_type>(
      do_not_init_tag("process symmetries"), process_symmetries.size(), max_symmetries + 1);
  auto process_symmetries_host_view = Kokkos::create_mirror_view(host_device_type(), process_symmetries_view);
  {
    for (ordinal_type i = 0, iend = process_symmetries_host_view.extent(0); i < iend; ++i) {
      const auto &symmetries = process_symmetries[i];
      process_symmetries_host_view(i, 0) = symmetries.size();
      for (ordinal_type j = 0, jend = symmetries.size(); j < jend; ++j)
        process_symmetries_host_view(i, j + 1) = symmetries[j];
    }
  }
  Kokkos::deep_copy(process_symmetries_view, process_symmetries_host_view);

  constraints_view = value_type_2d_view<ordinal_type, device_type>(
      do_not_init_tag("constraints"), processints_view.extent(0), variant_orderings_view.extent(0));
  auto constraints_host_view = Kokkos::create_mirror_view(host_device_type(), constraints_view);
  {
    using range_policy_type = Kokkos::RangePolicy<typename host_device_type::execution_space>;
    range_policy_type policy(0, processints_host_view.extent(0));
    Kokkos::parallel_for(policy, [&](const ordinal_type eid) {
      const ordinal_type conf_a_idx = processints_host_view(eid, 0);
      const ordinal_type conf_b_idx = processints_host_view(eid, 1);
      const ordinal_type process_idx = processints_host_view(eid, 2);
      const ordinal_type n_variants = variant_orderings_host_view.extent(0);
      const ordinal_type n_sites_conf = variant_orderings_host_view.extent(1);

      const auto &constraints = process_constraints[process_idx];

      using key_type = site_type[128];
      key_type a, a_variant, b, b_variant;
      for (ordinal_type k = 0; k < n_sites_conf; ++k) {
        a[k] = configurations_host_view(conf_a_idx, k);
        b[k] = configurations_host_view(conf_b_idx, k);
      }

      for (ordinal_type i = 0; i < n_variants; ++i) {
        constraints_host_view(eid, i) = -1; /// initialize with an invalid value
        for (ordinal_type k = 0; k < n_sites_conf; ++k)
          a_variant[k] = a[variant_orderings_host_view(i, k)];

        /// check if input variant meets the constraints
        bool input_conf_variant_valid(true);
        for (ordinal_type l = 0, lend = constraints.size(); l < lend; ++l) {
          const auto &constraint = constraints[l];
          if (a_variant[constraint[0]] != constraint[1]) {
            input_conf_variant_valid &= false;
          }
        }

        if (input_conf_variant_valid) {
          /// apply constraints to a_variant; b_variant matches to a_variant with constraints
          for (ordinal_type l = 0, lend = constraints.size(); l < lend; ++l) {
            const auto &constraint = constraints[l];
            a_variant[constraint[0]] = constraint[2];
          }

          ordinal_type n_variants_satisfying_constraints(0);
          for (ordinal_type j = 0; j < n_variants; ++j) {
            for (ordinal_type k = 0; k < n_sites_conf; ++k)
              b_variant[k] = b[variant_orderings_host_view(j, k)];

            /// check if the final conf variant meets the constraints
            bool constraints_satisfied(true);
            for (ordinal_type k = 0; k < n_sites_conf; ++k)
              constraints_satisfied &= (b_variant[k] == a_variant[k]);

            if (constraints_satisfied) {
              if (n_variants_satisfying_constraints == 0)
                constraints_host_view(eid, i) = j;
              ++n_variants_satisfying_constraints;
            }
          }
        }
      }
    });
  }
  if (verbose > 2) {
    std::cout << "\n";
    std::cout << "constraints : " << constraints_host_view.extent(0) << ", " << constraints_host_view.extent(1) << "\n";
    for (ordinal_type i = 0, iend = constraints_host_view.extent(0); i < iend; ++i) {
      ordinal_type sum = (constraints_host_view(i, 0) < 0);
      std::cout << "  instance: " << i << ", " << processes[processints_host_view(i, 2)] << ", ["
                << constraints_host_view(i, 0);
      for (ordinal_type j = 1, jend = constraints_host_view.extent(1); j < jend; ++j) {
        std::cout << ", " << constraints_host_view(i, j);
        sum += (constraints_host_view(i, j) < 0);
      }
      std::cout << "]\n";
      if (sum == ordinal_type(constraints_host_view.extent(1)))
        std::cout << " <------------- this event is null event\n";
    }
  }
  Kokkos::deep_copy(constraints_view, constraints_host_view);

  return 0;
}

ordinal_type validateOverrideInput(const boost::json::value &jv, boost::json::object &root,
                                   const ordinal_type verbose) {

  const FunctionScope func_scope("validateOverrideInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &root = jv.as_object();
    {
      const auto proc_it = root.find("processes");
      if (proc_it == root.end()) {
        // Do Nothing
      } else {
        const auto proc_rates_it = root.find("process rates");
        KINCAT_CHECK_ERROR(
            proc_rates_it == root.end(),
            "Error: [process rates] key does not exist in rates override file, even though processes are specified");
        const std::vector<std::string> proc_string = boost::json::value_to<std::vector<std::string>>(proc_it->value());
        KINCAT_CHECK_ERROR(proc_string.size() == 0, "Error: failed to find override process list");
        const auto proc_rates_root = proc_rates_it->value().as_object();
        ordinal_type rate_idx = 0;
        const std::vector<real_type> rates_data =
            boost::json::value_to<std::vector<real_type>>(proc_rates_root.find(std::to_string(rate_idx))->value());
        const auto n_samp_it = root.find("number of samples");
        if (n_samp_it == root.end()) {
          std::cout << "Failed to find the number of samples\n";
        }
        ordinal_type n_samp = boost::json::value_to<ordinal_type>(n_samp_it->value());
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught validating input\n" << e.what() << "\n";
  }
  return r_val;
}

ordinal_type validateRatesInput(const boost::json::value &jv, boost::json::object &root, const ordinal_type verbose) {

  const FunctionScope func_scope("validateRatesInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("rates");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }

    {
      const auto default_it = root.find("default rate");
      KINCAT_CHECK_ERROR(default_it == root.end(), "Error: [default rate] key does not exist in rates object");

      const auto proc_it = root.find("process specific rates");
      const auto instance_it = root.find("process instance specific rates");
      if (proc_it == root.end() && instance_it == root.end()) {
        std::cout << "Warning: Neither process nor instance specific rates exist in rates object \n All rates set to "
                     "default \n";
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught validating input\n" << e.what() << "\n";
  }

  return r_val;
}

ordinal_type parseRatesInput(const boost::json::object &root, const std::vector<std::string> &processes,
                             const value_type_2d_view<ordinal_type, device_type> &processints_view,
                             value_type_2d_view<real_type, device_type> &rates_view, const ordinal_type verbose) {
  const FunctionScope func_scope("parseRatesInput", __FILE__, __LINE__, verbose);

  /// parse data
  real_type default_rate;
  std::string type_string;
  std::vector<real_type> rates_data;

  /// view handling
  const auto processints_view_host = Kokkos::create_mirror_view_and_copy(host_device_type(), processints_view);
  rates_view = value_type_2d_view<real_type, device_type>(do_not_init_tag("rates_view"), 1, processints_view.extent(0));
  auto rates_view_host = Kokkos::create_mirror_view(host_device_type(), rates_view);

  try {
    //ordinal_type r_val(0);
    default_rate = boost::json::value_to<real_type>(root.find("default rate")->value());

    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "type : " << type_string << "\n";
      std::cout << "default rate : " << default_rate << "\n";
    }

    /// set the default value
    Kokkos::deep_copy(rates_view_host, default_rate);

    std::vector<real_type> rates_data(processes.size(), default_rate);

    const auto it = root.find("process specific rates");
    if (!(it == root.end())) {
      const auto rates_obj = root.find("process specific rates")->value().as_object();
      /// map and rates data
      std::map<std::string, std::vector<ordinal_type>> idxmap;
      for (ordinal_type i = 0, iend = processes.size(); i < iend; ++i) {
        idxmap[processes[i]].push_back(i);
      }

      std::vector<std::string> process_vector;
      std::vector<real_type> data_vector;
      for (const auto &key : processes) {
        const auto process_it = rates_obj.find(key);
        if (!(process_it == rates_obj.end())) {
          process_vector.push_back(key);
          const real_type temp_rate = boost::json::value_to<real_type>(process_it->value());
          data_vector.push_back(temp_rate);
        }
      }

      KINCAT_CHECK_ERROR(process_vector.size() != data_vector.size(),
                         "Error: process array length does not match to data array length");

      for (ordinal_type i = 0, iend = process_vector.size(); i < iend; ++i) {
        const auto &process_string = process_vector[i];
        const auto &idx = idxmap[process_string];
        for (ordinal_type j = 0, jend = idx.size(); j < jend; ++j)
          rates_data[idx[j]] = data_vector[i];
      }
      /// populate instance rates
      using range_policy_type = Kokkos::RangePolicy<typename host_device_type::execution_space>;
      range_policy_type policy(0, processints_view_host.extent(0));
      Kokkos::parallel_for(policy, [&](const ordinal_type i) {
        const ordinal_type idx = processints_view_host(i, 2);
        rates_view_host(0, i) = rates_data[idx];
      });
    }

    const auto instance_it = root.find("process instance specific rates");
    if (!(instance_it == root.end())) {
      const auto processints_obj = root.find("process instance specific rates")->value().as_object();

      // re-organize dictionary as coupled lists
      std::vector<ordinal_type> processints_vector;
      std::vector<real_type> data_vector;
      for (ordinal_type e_idx = 0, e_idxend = rates_view_host.extent(1); e_idx < e_idxend; e_idx++) {
        const auto e_it = processints_obj.find(std::to_string(e_idx));
        if (!(e_it == processints_obj.end())) {
          processints_vector.push_back(e_idx);
          const real_type temp_rate = boost::json::value_to<real_type>(e_it->value());
          data_vector.push_back(temp_rate);
        }
      }

      /// populate instance rates
      using range_policy_type = Kokkos::RangePolicy<typename host_device_type::execution_space>;
      range_policy_type policy(0, processints_vector.size());
      Kokkos::parallel_for(policy, [&](const ordinal_type i) {
        const ordinal_type idx = processints_vector[i];
        rates_view_host(0, idx) = data_vector[i];
      });
    }

    if (verbose > 2) {
      std::cout << "\n";
      std::cout << "rates : (" << rates_view_host.extent(0) << ", " << rates_view_host.extent(1) << ")\n";
      for (ordinal_type i = 0, iend = rates_view_host.extent(1); i < iend; ++i)
        std::cout << "[ " << i << ", " << rates_view_host(0, i) << "]\n";
    }

  } catch (const std::exception &e) {
    std::cout << "Execption " << __LINE__ << '\n';
    std::cerr << "Error: exception is caught parsing json input rates file\n" << e.what() << "\n";
  }

  Kokkos::deep_copy(rates_view, rates_view_host);

  return 0;
}

ordinal_type sortAndRemapDictionary(const value_type_2d_view<site_type, device_type> &configurations_view,
                                    const value_type_1d_view<ordinal_type, device_type> &configuration_sets_view,
                                    const value_type_2d_view<ordinal_type, device_type> &processints_view,
                                    const value_type_2d_view<ordinal_type, device_type> &constraints_view,
                                    const value_type_2d_view<real_type, device_type> &rates_view,
                                    const ordinal_type verbose) {
  auto configurations_host_view = Kokkos::create_mirror_view_and_copy(host_device_type(), configurations_view);
  auto configuration_sets_host_view = Kokkos::create_mirror_view_and_copy(host_device_type(), configuration_sets_view);
  std::vector<ordinal_type> idx_configurations, inverse_idx_configurations;
  std::vector<ordinal_type> idx_sets;
  std::vector<ordinal_type> idx_sets_configurations;
  {
    idx_configurations.resize(configurations_host_view.extent(0));
    std::iota(idx_configurations.begin(), idx_configurations.end(), 0);
    ordinal_type n_configs = configuration_sets_view.extent(0) - 1;
    
    for (ordinal_type seti = 0; seti < n_configs; seti++) {
      idx_sets.clear();
      for (ordinal_type i = ordinal_type(configuration_sets_host_view(seti)); i < ordinal_type(configuration_sets_host_view(seti+1)); i++) {
        idx_sets.push_back(idx_configurations[i]);
      }

      std::sort(idx_sets.begin(), idx_sets.end(),
              [&configurations_host_view](ordinal_type a, ordinal_type b) {
                for (ordinal_type i = 0, iend = configurations_host_view.extent(1); i < iend; ++i) {
                  if (configurations_host_view(a, i) == configurations_host_view(b, i))
                    continue;
                  return configurations_host_view(a, i) < configurations_host_view(b, i);
                }
                return false;
              });
      for (ordinal_type i = 0; i < idx_sets.size(); i++) {
        idx_sets_configurations.push_back(idx_sets[i]);
      }
    }
    permute(idx_sets_configurations, configurations_host_view);
    inverse_idx_configurations.resize(configurations_host_view.extent(0));
    for (ordinal_type i = 0, iend = inverse_idx_configurations.size(); i < iend; ++i) {
      inverse_idx_configurations[idx_sets_configurations[i]] = i;
    }
  }
  Kokkos::deep_copy(configurations_view, configurations_host_view);

  auto processints_host_view = Kokkos::create_mirror_view_and_copy(host_device_type(), processints_view);
  {
    for (ordinal_type i = 0, iend = processints_host_view.extent(0); i < iend; ++i) {
      const ordinal_type idx_conf_a = processints_host_view(i, 0);
      const ordinal_type idx_conf_b = processints_host_view(i, 1);
      processints_host_view(i, 0) = inverse_idx_configurations[idx_conf_a];
      processints_host_view(i, 1) = inverse_idx_configurations[idx_conf_b];
    }
  }

  std::vector<ordinal_type> idx_processints;
  {
    idx_processints.resize(processints_host_view.extent(0));
    std::iota(idx_processints.begin(), idx_processints.end(), 0);
    std::sort(idx_processints.begin(), idx_processints.end(), [&processints_host_view](ordinal_type a, ordinal_type b) {
      for (ordinal_type i = 0, iend = processints_host_view.extent(1); i < iend; ++i) {
        if (processints_host_view(a, i) == processints_host_view(b, i))
          continue;
        return processints_host_view(a, i) < processints_host_view(b, i);
      }
      return false;
    });
    permute(idx_processints, processints_host_view);
  }
  Kokkos::deep_copy(processints_view, processints_host_view);

  auto constraints_host_view = Kokkos::create_mirror_view_and_copy(host_device_type(), constraints_view);
  { permute(idx_processints, constraints_host_view); }
  Kokkos::deep_copy(constraints_view, constraints_host_view);

  auto rates_host_view = Kokkos::create_mirror_view_and_copy(host_device_type(), rates_view);
  {
    const value_type_1d_view<real_type, host_device_type> tmp(rates_host_view.data(), rates_host_view.extent(1));
    permute(idx_processints, tmp);
  }
  Kokkos::deep_copy(rates_view, rates_host_view);

  return 0;
}

ordinal_type validateSitesInput(const boost::json::value &jv, boost::json::object &root, const ordinal_type verbose) {

  const FunctionScope func_scope("validateSitesInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("sites");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }

    {
      const auto type_it = root.find("type");
      KINCAT_CHECK_ERROR(type_it == root.end(), "Error: [type] key does not exist in sites object");
      const std::string type_string = boost::json::value_to<std::string>(type_it->value());
      KINCAT_CHECK_ERROR(type_string != "random" && type_string != "file",
                         "Error: type string is not supported, available options are [random], [file]");

      std::vector<std::string> keys;
      keys.push_back("lattice");
      keys.push_back("random seed");
      keys.push_back("random fill ratio");
      keys.push_back("filename");

      for (const auto &key : keys) {
        const auto it = root.find(key);
        std::stringstream ss;
        ss << "Error: [" << key << "] key does not exist in sites object";
        KINCAT_CHECK_ERROR(it == root.end(), ss.str().c_str());
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught validating input\n" << e.what() << "\n";
  }

  return r_val;
}

ordinal_type parseSitesInput(const boost::json::object &root, std::string &site_init_type, ordinal_type &n_cells_x,
                             ordinal_type &n_cells_y, ordinal_type &n_basis, std::vector<real_type> &site_random_fill_ratio,
                             ordinal_type &site_random_seed, value_type_2d_view<site_type, device_type> &sites_view,
                             value_type_1d_view<real_type, device_type> &t_sites, const ordinal_type verbose) {
  const FunctionScope func_scope("parseSitesInput", __FILE__, __LINE__, verbose);
  /// parse data
  //real_type default_rate;
  std::string type_string;
  std::vector<real_type> rates_data;

  try {
    ordinal_type r_val(0);

    {
      std::vector<ordinal_type> tmp = boost::json::value_to<std::vector<ordinal_type>>(root.find("lattice")->value());
      n_cells_x = tmp[0];
      n_cells_y = tmp[1];
      n_basis = tmp[2];
    }

    site_init_type = boost::json::value_to<std::string>(root.find("type")->value());
    if (site_init_type == "random") {
      site_random_fill_ratio = boost::json::value_to<std::vector<real_type>>(root.find("random fill ratio")->value());
      site_random_seed = boost::json::value_to<ordinal_type>(root.find("random seed")->value());
    } else if (site_init_type == "file") {
      boost::json::object sites_root;
      {
        const std::string filename =
            replaceEnvironmentVariable(boost::json::value_to<std::string>(root.find("filename")->value()));
        boost::json::value jv;
        std::cout.flush();
        r_val = parseInputfile(filename, jv, verbose);
        const boost::json::object &tree = jv.as_object();
        const auto it = tree.find("sites");
        if (it == tree.end()) {
          sites_root = tree;
        } else {
          sites_root = it->value().as_object();
        }
      }
      //Check for number of samples in file
      ordinal_type n_samples_in = 0;
      ordinal_type batch_size_max = 100000; //This is intended to be far larger than any batch actually run. Will cause an error if matched. 
      for (ordinal_type i = 0, iend = batch_size_max; i < iend; i++) {
        const auto it = sites_root.find(std::to_string(i));
        if (it == sites_root.end()) {
          n_samples_in = i;
          break;
        } 
      }
      KINCAT_CHECK_ERROR((n_samples_in == 0), "WARNING:: Did not find samples in sites read-in file. Could be incorrect format.\n");
      t_sites = value_type_1d_view<real_type, device_type>(do_not_init_tag("t_sites"), n_samples_in);
      auto t_sites_host_view = Kokkos::create_mirror_view(host_device_type(), t_sites);
      sites_view =
          value_type_2d_view<site_type, device_type>(do_not_init_tag("sites"), 1, n_cells_x * n_cells_y * n_basis * n_samples_in);
      const auto sites_host_view = Kokkos::create_mirror_view(host_device_type(), sites_view);
      
      auto it = sites_root.find(std::to_string(0));
      boost::json::object sample_obj = it->value().as_object();
      real_type t_tmp = boost::json::value_to<real_type>(sample_obj.find("time")->value());
      std::vector<ordinal_type> tmp =
        boost::json::value_to<std::vector<ordinal_type>>(sample_obj.find("data")->value());
      KINCAT_CHECK_ERROR((tmp.size() != (n_cells_x * n_cells_y * n_basis)), "Error: sites file sample does not match to lattice information");
      std::vector<ordinal_type> complete_tmp = tmp;
      std::vector<real_type> t_sites_vector ;
      t_sites_vector.push_back(t_tmp);

      for (ordinal_type i = 1; i < n_samples_in; i++) {
        const auto it = sites_root.find(std::to_string(i));
        sample_obj = it->value().as_object();
        real_type t_tmp = boost::json::value_to<real_type>(sample_obj.find("time")->value());
        t_sites_vector.push_back(t_tmp);
        std::vector<ordinal_type> tmp =
          boost::json::value_to<std::vector<ordinal_type>>(sample_obj.find("data")->value());
        KINCAT_CHECK_ERROR((tmp.size() != (n_cells_x * n_cells_y * n_basis)), "Error: sites file sample does not match to lattice information");
        complete_tmp.insert(complete_tmp.end(), tmp.begin(), tmp.end());
      }

      KINCAT_CHECK_ERROR(sites_host_view.extent(1) != complete_tmp.size(),
                         "Error: sites file samples do not match to lattice information");
      copyStdVectorToKokkosView(t_sites_vector, t_sites_host_view);
      Kokkos::deep_copy(t_sites, t_sites_host_view);
      Kokkos::deep_copy(Kokkos::subview(sites_host_view, 0, Kokkos::ALL()),
                        value_type_1d_view<ordinal_type, host_device_type>(complete_tmp.data(), complete_tmp.size()));
      Kokkos::deep_copy(sites_view, sites_host_view);
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: exception is caught parsing json input file sites object\n" << e.what() << "\n";
  }

  return 0;
}

ordinal_type validateDumpInput(const boost::json::value &jv, boost::json::object &root, const ordinal_type verbose) {

  const FunctionScope func_scope("validateDumpInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("dump");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }
    {
      const auto file_it = root.find("dump filename");
      KINCAT_CHECK_ERROR(file_it == root.end(), "Error: [dump filename] key does not exist in dump object");
      const std::string file_string = boost::json::value_to<std::string>(file_it->value());
      const auto int_it = root.find("dump interval");
      if (!(int_it == root.end())) { // output interval included
        const real_type dump_interval = boost::json::value_to<real_type>(int_it->value());
        KINCAT_CHECK_ERROR((dump_interval < 0), "Error: [dump interval] may not be negative")
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught parsing json input statistics section\n" << e.what() << "\n";
  }
  return r_val;
}

ordinal_type validateStatsInput(const boost::json::value &jv, boost::json::object &root, const ordinal_type verbose) {

  const FunctionScope func_scope("validateStatsInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("statistics");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }

    {

      const auto file_it = root.find("stats filename");
      KINCAT_CHECK_ERROR(file_it == root.end(), "Error: [output filename] key does not exist in statistics object");
      const std::string file_string = boost::json::value_to<std::string>(file_it->value());
      const auto types_it = root.find("types");
      KINCAT_CHECK_ERROR(types_it == root.end(), "Error: [types] key does not exist in statistics object");
      const std::vector<std::string> types_list = boost::json::value_to<std::vector<std::string>>(types_it->value());
      const auto int_it = root.find("stats interval");
      if (!(int_it == root.end())) { // output interval included
        const real_type stats_interval = boost::json::value_to<real_type>(int_it->value());
        KINCAT_CHECK_ERROR((stats_interval < 0), "Error: [output interval] may not be negative")
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught parsing json input statistics section\n" << e.what() << "\n";
  }
  return r_val;
}

ordinal_type validateSolverInput(const boost::json::value &jv, boost::json::object &root, const ordinal_type verbose) {

  const FunctionScope func_scope("validateSolverInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("solver");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }

    {
      const auto type_it = root.find("type");
      KINCAT_CHECK_ERROR(type_it == root.end(), "Error: [type] key does not exist in solver object");
      const std::string type_string = boost::json::value_to<std::string>(type_it->value());
      KINCAT_CHECK_ERROR(!(type_string == "serial-rta" || type_string == "sublattice" || type_string == "batch-rta" ||
                           type_string == "timewarp"),
                         "Error: type string is not supported, available options are [serial-rta], [sublattice]");
      if (type_string == "sublattice" || type_string == "timewarp") {
        const auto it = root.find("domain");
        std::stringstream ss;
        ss << "Error: [domain] key needed in solver object for parallel types";
        KINCAT_CHECK_ERROR(it == root.end(), ss.str().c_str());
      }

      std::vector<std::string> keys;
      keys.push_back("random seed");
      keys.push_back("random pool size");
      keys.push_back("max number of kmc kernel launches");
      keys.push_back("max number of kmc steps per kernel launch");
      keys.push_back("time range");

      for (const auto &key : keys) {
        const auto it = root.find(key);
        std::stringstream ss;
        ss << "Error: [" << key << "] key does not exist in solver object";
        KINCAT_CHECK_ERROR(it == root.end(), ss.str().c_str());
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught validating input\n" << e.what() << "\n";
  }

  return r_val;
}

ordinal_type parseSolverInput(const boost::json::object &root, std::string &solver_type, ordinal_type &n_cells_domain_x,
                              ordinal_type &n_cells_domain_y, ordinal_type &solver_random_seed,
                              ordinal_type &solver_random_pool_size, ordinal_type &n_kmc_kernel_launches,
                              ordinal_type &n_kmc_steps_per_kernel_launch, real_type &t_begin, real_type &t_end,
                              real_type &t_dt, const ordinal_type verbose) {
  const FunctionScope func_scope("parseSolverInput", __FILE__, __LINE__, verbose);

  /// parse data
  real_type default_rate;
  std::string type_string;
  std::vector<real_type> rates_data;

  try {
    ordinal_type r_val(0);
    solver_type = boost::json::value_to<std::string>(root.find("type")->value());
    std::cout << "solver_type = " << solver_type << '\n';
    if (solver_type == "sublattice" || solver_type == "timewarp") {
      std::vector<ordinal_type> tmp = boost::json::value_to<std::vector<ordinal_type>>(root.find("domain")->value());
      n_cells_domain_x = tmp[0];
      n_cells_domain_y = tmp[1];
    }
    
    solver_random_seed = boost::json::value_to<ordinal_type>(root.find("random seed")->value());
    solver_random_pool_size = boost::json::value_to<ordinal_type>(root.find("random pool size")->value());
    n_kmc_kernel_launches =
        boost::json::value_to<ordinal_type>(root.find("max number of kmc kernel launches")->value());
    n_kmc_steps_per_kernel_launch =
        boost::json::value_to<ordinal_type>(root.find("max number of kmc steps per kernel launch")->value());
    {
      std::vector<real_type> tmp = boost::json::value_to<std::vector<real_type>>(root.find("time range")->value());
      t_begin = tmp[0];
      t_end = tmp[1];
      if (tmp.size() < 3)
        t_dt = 0;
      else
        t_dt = tmp[2];
    }

  } catch (const std::exception &e) {
    std::cout << "Execption " << __LINE__ << '\n';
    std::cerr << "Error: exception is caught parsing json input solver object\n" << e.what() << "\n";
  }

  return 0;
}
ordinal_type validateEnsembleInput(const boost::json::value &jv, boost::json::object &root,
                                   const ordinal_type verbose) {

  const FunctionScope func_scope("validateEnsembleInput", __FILE__, __LINE__, verbose);

  ordinal_type r_val(0);

  try {
    const boost::json::object &tree = jv.as_object();
    {
      const auto it = tree.find("ensemble");
      if (it == tree.end()) {
        root = tree;
      } else {
        root = it->value().as_object();
      }
    }

    {
      {
        /// required keys
        std::vector<std::string> keys;
        keys.push_back("number of samples");
        for (const auto &key : keys) {
          const auto it = root.find(key);
          std::stringstream ss;
          ss << "Error: [" << key << "] key does not exist in solver object";
          KINCAT_CHECK_ERROR(it == root.end(), ss.str().c_str());
        }
      }
      {
        /// TODO:: test optional keys
      }
    }
  } catch (const std::exception &e) {
    r_val = -1;
    std::cerr << "Error: exception is caught validating input\n" << e.what() << "\n";
  }

  return r_val;
}

ordinal_type parseEnsembleInput(const boost::json::object &root, ordinal_type &n_samples,
                                bool &is_solver_random_number_variation, std::vector<real_type> &site_random_fill_ratio,
                                std::vector<std::string> &process_override,
                                std::vector<real_type> &process_rates_override,
                                std::vector<ordinal_type> &instance_override,
                                std::vector<real_type> &instance_rates_override, const ordinal_type verbose) {
  const FunctionScope func_scope("parseEnsembleInput", __FILE__, __LINE__, verbose);

  /// parse data
  try {
    ordinal_type r_val(0);
    n_samples = boost::json::value_to<ordinal_type>(root.find("number of samples")->value());
    {
      const auto it = root.find("solver random number variations");
      if (it == root.end()) {
        if (verbose) {
          std::cout << "key: solver random number variations does not exist\n";
        }
      } else {
        boost::json::object solver_random_number_variation_root = it->value().as_object();
        const std::string apply =
            boost::json::value_to<std::string>(solver_random_number_variation_root.find("apply")->value());
        is_solver_random_number_variation = (apply == "enabled");
      }
    }
    {
      const auto it = root.find("sites random variations");
      if (it == root.end()) {
        if (verbose) {
          std::cout << "key: sites random variations does not exist\n";
        }
      } else {
        boost::json::object sites_random_variation_root = it->value().as_object();
        const std::string apply =
            boost::json::value_to<std::string>(sites_random_variation_root.find("apply")->value());
        site_random_fill_ratio.clear();
        if (apply == "enabled") {
          std::cout << __LINE__ << "\n";
          const auto fill_it = sites_random_variation_root.find("random fill ratio");
          if (fill_it != sites_random_variation_root.end()) { // fill ratio variants included
            std::vector<real_type> tmp;
            boost::json::object fill_root = fill_it->value().as_object();
            ordinal_type tmp_size = 0;
            for (ordinal_type fill_idx = 0; fill_idx < n_samples; fill_idx++) {
              KINCAT_CHECK_ERROR((fill_root.find(std::to_string(fill_idx)) == fill_root.end()), 
                                  "Missing fill ratios sample in 'random fill ratio'");
              tmp = boost::json::value_to<std::vector<real_type>>(fill_root.find(std::to_string(fill_idx))->value());
              if (tmp_size == 0) {
                tmp_size = tmp.size();
              }
              KINCAT_CHECK_ERROR((tmp.size() != (tmp_size)),
                                  "Sample fill ratio does not have same number of entries (species) as before");
              for (ordinal_type tmp_idx = 0; tmp_idx < tmp.size(); tmp_idx++) {
                site_random_fill_ratio.push_back(tmp[tmp_idx]);
              }
            }
          } else {// each sample will be filled by the ratios given in the sites object.
          }
        }
      }
    }
    {
      const auto it = root.find("rates variations");
      if (it == root.end()) {
        if (verbose) {
          std::cout << "key: rates variations does not exist\n";
        }
      } else {
        boost::json::object rates_variation_root = it->value().as_object();
        const std::string apply = boost::json::value_to<std::string>(rates_variation_root.find("apply")->value());
        process_override.clear();
        process_rates_override.clear();
        instance_override.clear();
        instance_rates_override.clear();
        if (apply == "enabled") {
          const std::string type = boost::json::value_to<std::string>(rates_variation_root.find("type")->value());
          if (type == "inlined") {
            const auto process_rates_it = rates_variation_root.find("processes");
            if (process_rates_it == rates_variation_root.end()) {
              if (verbose) {
                std::cout << "key: processes in rates variations does not exist\n";
              }
            } else {
              process_override =
                  boost::json::value_to<std::vector<std::string>>(rates_variation_root.find("processes")->value());
              ordinal_type n_proc_override = process_override.size();
              const auto rates_it = rates_variation_root.find("process rates");
              KINCAT_CHECK_ERROR((rates_it == rates_variation_root.end()),
                                 "no process rates given, even though processes given in 'rates variations' object");
              std::vector<real_type> rates_vector, row_vector;
              boost::json::object rates_root = rates_it->value().as_object();
              for (ordinal_type rate_idx = 0, rate_idx_end = n_samples; rate_idx < rate_idx_end; rate_idx++) {
                KINCAT_CHECK_ERROR((rates_root.find(std::to_string(rate_idx)) == rates_root.end()),
                                   "Missing process rates sample in 'rates variations'");
                row_vector =
                    boost::json::value_to<std::vector<real_type>>(rates_root.find(std::to_string(rate_idx))->value());
                KINCAT_CHECK_ERROR((row_vector.size() != n_proc_override),
                                   "Incorrect number of rates per sample in 'rates variations' ");
                for (ordinal_type proc_idx = 0, proc_idx_end = n_proc_override; proc_idx < proc_idx_end; proc_idx++) {
                  rates_vector.push_back(row_vector[proc_idx]);
                }
              }
              process_rates_override = rates_vector;
            }
            const auto instance_rates_it = rates_variation_root.find("process instances");
            if (instance_rates_it == rates_variation_root.end()) {
              if (verbose) {
                std::cout << "key: process instances in rates variations does not exist\n";
              }
            } else {
              // Clunky but works, dictionary inputs
              instance_override = boost::json::value_to<std::vector<ordinal_type>>(
                  rates_variation_root.find("process instances")->value());
              ordinal_type n_inst_override = instance_override.size();
              const auto inst_it = rates_variation_root.find("instance rates");
              KINCAT_CHECK_ERROR((inst_it == rates_variation_root.end()),
                                 "no instance rates given, even though instances given in 'rates variations' object");
              std::vector<real_type> inst_rates_vector, row_vector;
              boost::json::object inst_root = inst_it->value().as_object();
              for (ordinal_type rate_idx = 0, rate_idx_end = n_samples; rate_idx < rate_idx_end; rate_idx++) {
                KINCAT_CHECK_ERROR((inst_root.find(std::to_string(rate_idx)) == inst_root.end()),
                                   "Missing instance sample in 'rates variations'");
                row_vector =
                    boost::json::value_to<std::vector<real_type>>(inst_root.find(std::to_string(rate_idx))->value());
                KINCAT_CHECK_ERROR((row_vector.size() != n_inst_override),
                                   "Incorrect number of instance rates per sample in 'rates variations' ");
                for (ordinal_type inst_idx = 0, inst_idx_end = n_inst_override; inst_idx < inst_idx_end; inst_idx++) {
                  inst_rates_vector.push_back(row_vector[inst_idx]);
                }
              }
              instance_rates_override = inst_rates_vector;
            }
            KINCAT_CHECK_ERROR((process_rates_override.size() == 0 && instance_rates_override.size() == 0),
                               "Rates variation enabled, but no rate variations given.")

          } else if (type == "file") {
            std::string override_file =
                boost::json::value_to<std::string>(rates_variation_root.find("override filename")->value());
            if (verbose > 0) {
              std::cout << "override_file: " << override_file << '\n';
            }
            boost::json::value jv;
            r_val = KinCat::parseInputfile(override_file, jv, verbose);
            KINCAT_CHECK_ERROR(r_val, "Error: fails to parse override rates file");

            boost::json::object override_root;
            r_val = KinCat::validateOverrideInput(jv, override_root, verbose);
            KINCAT_CHECK_ERROR(r_val, "Error: fails to validate rate override file");

            override_root = jv.as_object();
            if (override_root.find("number of samples") != override_root.end()) {
              ordinal_type n_samples_check =
                  boost::json::value_to<ordinal_type>(override_root.find("number of samples")->value());
              if (n_samples_check != n_samples) {
                std::cout << "WARNING: number of samples specified in rates variations override file does not match "
                             "the number specified in ensemble object!!! \n";
              }
            }
            const auto process_rates_it = override_root.find("processes");
            if (process_rates_it == override_root.end()) {
              if (verbose) {
                std::cout << "key: processes in rates variations does not exist\n";
              }
            } else {
              process_override =
                  boost::json::value_to<std::vector<std::string>>(override_root.find("processes")->value());
              ordinal_type n_proc_override = process_override.size();
              std::cout << "n_proc_override: " << n_proc_override << '\n';
              const auto rates_it = override_root.find("process rates");
              KINCAT_CHECK_ERROR((rates_it == override_root.end()),
                                 "no process rates given, even though processes given in 'rates variations' object");
              std::vector<real_type> rates_vector, row_vector;
              boost::json::object rates_root = rates_it->value().as_object();
              for (ordinal_type rate_idx = 0, rate_idx_end = n_samples; rate_idx < rate_idx_end; rate_idx++) {
                KINCAT_CHECK_ERROR((rates_root.find(std::to_string(rate_idx)) == rates_root.end()),
                                   "Missing process rates sample in 'rates variations'");
                row_vector =
                    boost::json::value_to<std::vector<real_type>>(rates_root.find(std::to_string(rate_idx))->value());
                KINCAT_CHECK_ERROR((row_vector.size() != n_proc_override),
                                   "Incorrect number of rates per sample in 'rates variations' ");
                for (ordinal_type proc_idx = 0, proc_idx_end = n_proc_override; proc_idx < proc_idx_end; proc_idx++) {
                  rates_vector.push_back(row_vector[proc_idx]);
                }
              }
              process_rates_override = rates_vector;
            }
            const auto instance_rates_it = override_root.find("process instances");
            if (instance_rates_it == override_root.end()) {
              if (verbose) {
                std::cout << "key: process instances in rates variations does not exist\n";
              }
            } else {
              // Clunky but works, dictionary inputs
              instance_override =
                  boost::json::value_to<std::vector<ordinal_type>>(override_root.find("process instances")->value());
              ordinal_type n_inst_override = instance_override.size();
              const auto inst_it = override_root.find("instance rates");
              KINCAT_CHECK_ERROR((inst_it == override_root.end()),
                                 "no instance rates given, even though instances given in 'rates variations' object");
              std::vector<real_type> inst_rates_vector, row_vector;
              boost::json::object inst_root = inst_it->value().as_object();
              for (ordinal_type rate_idx = 0, rate_idx_end = n_samples; rate_idx < rate_idx_end; rate_idx++) {
                KINCAT_CHECK_ERROR((inst_root.find(std::to_string(rate_idx)) == inst_root.end()),
                                   "Missing instance sample in 'rates variations'");
                row_vector =
                    boost::json::value_to<std::vector<real_type>>(inst_root.find(std::to_string(rate_idx))->value());
                KINCAT_CHECK_ERROR((row_vector.size() != n_inst_override),
                                   "Incorrect number of instance rates per sample in 'rates variations' ");
                for (ordinal_type inst_idx = 0, inst_idx_end = n_inst_override; inst_idx < inst_idx_end; inst_idx++) {
                  inst_rates_vector.push_back(row_vector[inst_idx]);
                }
              }
              instance_rates_override = inst_rates_vector;
            }
            KINCAT_CHECK_ERROR((process_rates_override.size() == 0 && instance_rates_override.size() == 0),
                               "Rates variation enabled, but no rate variations given in rates override file.")
          }
        }
      }
    }

  } catch (const std::exception &e) {
    std::cout << "Execption " << __LINE__ << '\n';
    std::cerr << "Error: exception is caught parsing json input rates file\n" << e.what() << "\n";
  }

  return 0;
}

ordinal_type Input::parse(const std::string input_filename, const ordinal_type verbose) {
  ordinal_type r_val(0);

  try {
    /// top-level input parsing
    boost::json::value jv_input;
    r_val = KinCat::parseInputfile(input_filename, jv_input, verbose);
    KINCAT_CHECK_ERROR(r_val, "Error: fails to parse input json file [input]");

    boost::json::object jv_input_obj = jv_input.as_object();
    {
      const auto it = jv_input_obj.find("kincat");
      KINCAT_CHECK_ERROR(it == jv_input_obj.end(), "Error: fails to find key [kincat] from input json file");
      jv_input_obj = it->value().as_object();
    }

    /// dictionary
    {
      std::string filename;
      {
        const auto it = jv_input_obj.find("dictionary");
        KINCAT_CHECK_ERROR(it == jv_input_obj.end(), "Error: fails to find key [dictionary] in input object");
        filename = KinCat::replaceEnvironmentVariable(boost::json::value_to<std::string>(it->value()));
      }
      if (verbose > 0) {
        std::cout << "dictionary file : " << filename << "\n";
      }

      {
        boost::json::value jv;
        r_val = KinCat::parseInputfile(filename, jv, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to parse input json file [dictionary]");
        boost::json::object root;
        r_val = KinCat::validateDictionaryInput(jv, root, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to validate dictionary input");
        r_val = KinCat::parseDictionaryInput(
            root, _edge_vectors, _basis_vectors, _symmetry_operations, _site_coordinates, _n_species,
            _variant_orderings, _symmetry_orderings, _n_cells_interaction_x, _n_cells_interaction_y, _configurations, _configuration_sets, _processints,
            _constraints, _process_symmetries, _processes, _process_index, _process_constraints, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to parse dictionary input");
      }
    }

    /// rates
    {
      std::string filename;
      {
        const auto it = jv_input_obj.find("rates");
        KINCAT_CHECK_ERROR(it == jv_input_obj.end(), "Error: fails to find key [rates] in input object");
        filename = KinCat::replaceEnvironmentVariable(boost::json::value_to<std::string>(it->value()));
      }
      if (verbose > 0) {
        std::cout << "rates file : " << filename << "\n";
      }

      {
        boost::json::value jv;
        r_val = KinCat::parseInputfile(filename, jv, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to parse input json file [rates]");
        boost::json::object root;
        r_val = KinCat::validateRatesInput(jv, root, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to validate rates input");

        r_val = KinCat::parseRatesInput(root, _processes, _processints, _rates, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to parse rates input");
      }
    }

    /// sort and remap
    { r_val = KinCat::sortAndRemapDictionary(_configurations, _configuration_sets, _processints, _constraints, _rates); }

    /// sites
    {
      boost::json::value jv;
      {
        const auto it = jv_input_obj.find("sites");
        KINCAT_CHECK_ERROR(it == jv_input_obj.end(), "Error: fails to find key [sites] in input object");
        jv = it->value();
      }

      boost::json::object root;
      r_val = KinCat::validateSitesInput(jv, root, verbose);
      KINCAT_CHECK_ERROR(r_val, "Error: fails to validate sites input");
      r_val = KinCat::parseSitesInput(root,
                                      /// site init type
                                      _site_init_type,
                                      /// lattice
                                      _n_cells_x, _n_cells_y, _n_basis,
                                      /// random
                                      _site_random_fill_ratio, _site_random_seed,
                                      /// file
                                      _sites, _t_sites, verbose);
      KINCAT_CHECK_ERROR(_n_basis != _basis_vectors.extent(0),
                         "Error: number of basis does not match to basis vectors size");
    }

    /// dump file
    {
      {
        boost::json::value jv;
        const auto it = jv_input_obj.find("dump");
        if (it == jv_input_obj.end()) { // Dump option not included, ignore functionality
          std::string _dump = "none";
          _dump_filename = _dump;
        } else {
          jv = it->value();
          boost::json::object root;
          r_val = KinCat::validateDumpInput(jv, root, verbose);
          KINCAT_CHECK_ERROR(r_val, "fails to validate dump input object");
          const auto file_it = root.find("dump filename");
          _dump_filename = KinCat::replaceEnvironmentVariable(boost::json::value_to<std::string>(file_it->value()));
          const auto file_rs_it = root.find("restart filename");
          if (file_rs_it == root.end()) { //Restart file not specified
            std::string _save = "restart_sites.json";
            _restart_save_filename = _save; 
          } else {
            _restart_save_filename = KinCat::replaceEnvironmentVariable(boost::json::value_to<std::string>(file_rs_it->value()));
          }
          const auto int_it = root.find("dump interval");
          if (int_it == root.end()) {
            _dump_interval = 0.0;
          } else {
            _dump_interval = boost::json::value_to<real_type>(int_it->value());
            if (_dump_interval < 0.0) {
              _dump_interval = 0.0;
            }
          }
        }
      }
      if (verbose > 0) {
        std::cout << "dump file : " << _dump_filename << "\n";
        std::cout << "dump interval : " << _dump_interval << "\n";
      }
    }
    /// stats file
    {
      {
        boost::json::value jv;
        const auto it = jv_input_obj.find("statistics");
        if (it == jv_input_obj.end()) { // Statistics option not included, ignore functionality
          std::string _sts = "none";
          _stats_filename = _sts;
        }

        else {
          jv = it->value();
          boost::json::object root;
          r_val = KinCat::validateStatsInput(jv, root, verbose);
          KINCAT_CHECK_ERROR(r_val, "fails to validate statistics input object");
          const auto file_it = root.find("stats filename");
          _stats_filename = KinCat::replaceEnvironmentVariable(boost::json::value_to<std::string>(file_it->value()));
          const auto types_it = root.find("types");
          std::vector<std::string> types_list = boost::json::value_to<std::vector<std::string>>(types_it->value());
          _stats_list.clear();
          _stats_list = (boost::json::value_to<std::vector<std::string>>(types_it->value()));
          const auto int_it = root.find("stats interval");
          if (int_it == root.end()) {
            _stats_interval = 0.0;
          } else {
            _stats_interval = boost::json::value_to<real_type>(int_it->value());
            if (_stats_interval < 0.0) {
              _stats_interval = 0.0;
            }
          }
        }
      }
      if (verbose > 0) {
        std::cout << "statistics file : " << _stats_filename << "\n";
        std::cout << "statistics types: ";
        for (ordinal_type i = 0, iend = _stats_list.size(); i < iend; i++) {
          std::cout << ' ' << _stats_list[i];
        }
        std::cout << '\n';
      }
    }
    /// solver
    {
      boost::json::value jv;
      {
        const auto it = jv_input_obj.find("solver");
        KINCAT_CHECK_ERROR(it == jv_input_obj.end(), "Error: fails to find key [solver] in input object");
        jv = it->value();
      }

      boost::json::object root;
      r_val = KinCat::validateSolverInput(jv, root, verbose);
      KINCAT_CHECK_ERROR(r_val, "Error: fails to validate solver input");
      r_val = KinCat::parseSolverInput(root,
                                       /// solver type
                                       _solver_type,
                                       /// domain size per worker
                                       _n_cells_domain_x, _n_cells_domain_y,
                                       /// random variable control
                                       _solver_random_seed, _solver_random_pool_size,
                                       /// kmc steps
                                       _n_kmc_kernel_launches, _n_kmc_steps_per_kernel_launch,
                                       /// time range
                                       _t_begin, _t_end, _t_dt, verbose);
      if (_solver_type == "serial-rta" ||
          _solver_type == "batch-rta") { // Enforce the domain size being the same as the lattice size. This also
                                         // instantiates _n_cells_domain_x/y for serial solvers.
        _n_cells_domain_x = _n_cells_x;
        _n_cells_domain_y = _n_cells_y;
      }
      KINCAT_CHECK_ERROR(r_val, "Error: fails to validate solver input");
      if (_solver_type == "sublattice") {
        KINCAT_CHECK_ERROR((_n_cells_x / _n_cells_domain_x) % 2 != 0,
                           "Error: lattice is not an even multiple of domain x");
        KINCAT_CHECK_ERROR((_n_cells_y / _n_cells_domain_y) % 2 != 0,
                           "Error: lattice is not an even multiple of domain y");
      }
    }
    /// ensemble
    {
      const auto it = jv_input_obj.find("ensemble");
      if (it == jv_input_obj.end()) {
        /// do nothing, ensemble input is optional
      } else {
        boost::json::value jv = it->value();
        boost::json::object root;
        r_val = KinCat::validateEnsembleInput(jv, root, verbose);
        KINCAT_CHECK_ERROR(r_val, "Error: fails to validate ensemble input");
        r_val = KinCat::parseEnsembleInput(root, _n_samples, _is_solver_random_number_variation,
                                           _site_random_fill_ratio_array, _process_override_array,
                                           _process_rates_override_array, _instance_override_array,
                                           _instance_rates_override_array, verbose);

        if ((_n_samples > 1) && _site_random_fill_ratio_array.size() == 0 && _site_random_fill_ratio.size() != 0) { // Vary sites, but not fill ratios, sites not read-in
          ordinal_type n_species = _n_species - 1;
          _site_random_fill_ratio_array.resize(_n_samples * n_species);
          for (ordinal_type samp_idx = 0; samp_idx < _n_samples; samp_idx++) {
            for (ordinal_type spec_idx = 0; spec_idx < n_species; spec_idx++) {
              _site_random_fill_ratio_array[(samp_idx * n_species) + spec_idx] = _site_random_fill_ratio[spec_idx];
            }
          }
        }
        KINCAT_CHECK_ERROR(r_val, "Error: fails to validate solver input");
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: exception is caught parsing dictionary input\n" << e.what() << "\n";
    return -1;
  }
  return r_val;
}

} // namespace KinCat