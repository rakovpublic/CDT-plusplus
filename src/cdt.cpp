/*******************************************************************************
 Causal Dynamical Triangulations in C++ using CGAL

 Copyright © 2013–2026 Adam Getchell
 ******************************************************************************/

/// @file cdt.cpp
/// @brief The main executable
/// @author Adam Getchell
/// @details A program that generates spacetime ensembles. Inspired by
/// https://github.com/ucdavis/CDT.

#include <CGAL/Real_timer.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <boost/program_options.hpp>
#if defined(CDT_ENABLE_PARALLEL_TRIANGULATION) && \
    CDT_ENABLE_PARALLEL_TRIANGULATION
#include <oneapi/tbb/global_control.h>
#endif

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <Metropolis.hpp>
#include <optional>
#include <Simulation_output.hpp>
#include <string>
#include <string_view>
#include <utility>

#include "Runtime_config.hpp"
#include "Version.hpp"

using Timer = CGAL::Real_timer;

using namespace cdt;
using namespace std;
namespace po = boost::program_options;

/// Help text used by Boost.Program_options
static constexpr string_view USAGE{
    R"(Causal Dynamical Triangulations in C++ using CGAL.

Copyright (c) 2013-2026 Adam Getchell

A program that generates d-dimensional triangulated spacetimes
with a defined causal structure and evolves them according
to the Metropolis algorithm. Specify the number of passes to control
how much evolution is desired. Each pass attempts a number of ergodic
moves equal to the number of simplices in the simulation.

Usage: ./cdt ((--spherical | --toroidal) -n SIMPLICES -t TIMESLICES
             [-d DIM] [--init INITIAL RADIUS] [--foliate FOLIATION SPACING]
             | --input INITIAL.off)
            (-k K --alpha ALPHA --lambda LAMBDA
             | --kappa0 KAPPA0 --kappa4 KAPPA4 --Delta DELTA
               [--target-n4 TARGET] [--volume-epsilon EPSILON]
               [--thermalization STEPS] [--measurement-interval STEPS]
               [--chain-id ID] [--run-id ID] [--output-dir DIR])
            [--no-output] [--seed SEED] [--threads THREADS]
            [-p PASSES] [-c CHECKPOINT]
      ./cdt --resume CHECKPOINT.off [-p TOTAL_PASSES] [--no-output]

Optional arguments are in square brackets.

Examples:
./cdt --spherical -n 32000 -t 11 --alpha 0.6 -k 1.1 --lambda 0.1 --passes 1000
./cdt -s -n32000 -t11 -a.6 -k1.1 -l.1 -p1000 --seed 92
./cdt --input <initialize-output>.off -a.6 -k1.1 -l.1 -p1000 --seed 93
./cdt --resume <checkpoint>.off
./cdt --spherical -d4 -n32000 -t16 --kappa0 2.2 --kappa4 0.9 --Delta 0.6

Options)"};

/// @brief The main path of the CDT++ program
/// @param argc Argument count = 1 + number of arguments
/// @param argv Argument vector passed to Boost.Program_options
/// @return Integer value 0 if successful, 1 on failure
auto main(int const argc, char* const argv[]) -> int
try
{
  std::string const intro{USAGE};
  // Parsed arguments
  long long               simplices{};
  long long               timeslices{};
  long long               dimensions{};
  double                  initial_radius{};
  double                  foliation_spacing{};
  long double             alpha{};
  long double             k{};
  long double             lambda{};
  long double             kappa_0{};
  long double             kappa_4{};
  long double             Delta{};
  long double             volume_epsilon{};
  long long               passes{};
  long long               checkpoint{};
  std::uint64_t           seed{};
  long long               threads{};
  long long               target_N4{};
  long long               thermalization{};
  long long               measurement_interval{};
  std::string             input_path;
  std::string             resume_path;
  std::string             chain_id;
  std::string             run_id;
  std::string             output_dir;

  po::options_description description(intro);
  description.add_options()("help,h", "Show this message")(
      "version,v", "Show program version")("spherical,s", "Spherical topology")(
      "toroidal,e", "Toroidal topology")("simplices,n",
                                         po::value<long long>(&simplices),
                                         "Approximate number of simplices")(
      "timeslices,t", po::value<long long>(&timeslices),
      "Number of timeslices")(
      "dimensions,d", po::value<long long>(&dimensions)->default_value(3),
      "Dimensionality")("init,i",
                        po::value<double>(&initial_radius)->default_value(1.0),
                        "Initial radius")(
      "foliate,f", po::value<double>(&foliation_spacing)->default_value(1.0),
      "Foliation spacing")("input", po::value<std::string>(&input_path),
                           "Initial-triangulation payload (requires .meta)")(
      "resume", po::value<std::string>(&resume_path),
      "Resume the identical Markov chain from a checkpoint")(
      "no-output", "Do not write checkpoint or final triangulation files")(
      "seed", po::value<std::uint64_t>(&seed),
      "Root random seed (default: operating-system entropy)")(
      "threads", po::value<long long>(&threads)->default_value(1),
      "Maximum worker threads for supported Delaunay operations")(
      "alpha,a", po::value<long double>(&alpha),
      "Negative squared geodesic length of 1-d timelike edges")(
      "k,k", po::value<long double>(&k), "K = 1/(8*pi*G_newton)")(
      "lambda,l", po::value<long double>(&lambda), "K * Cosmological constant")(
      "kappa0", po::value<long double>(&kappa_0),
      "4D bare inverse Newton coupling")(
      "kappa4", po::value<long double>(&kappa_4),
      "4D bare cosmological coupling")(
      "Delta", po::value<long double>(&Delta), "4D asymmetry coupling")(
      "target-n4", po::value<long long>(&target_N4)->default_value(0),
      "4D fixed-volume target")(
      "volume-epsilon",
      po::value<long double>(&volume_epsilon)->default_value(0.0L),
      "Quadratic fixed-volume strength")(
      "thermalization",
      po::value<long long>(&thermalization)->default_value(0),
      "4D thermalization steps discarded before measurements")(
      "measurement-interval",
      po::value<long long>(&measurement_interval)->default_value(1),
      "4D measurement interval")(
      "chain-id", po::value<std::string>(&chain_id)->default_value("chain-0"),
      "Independent chain identifier")(
      "run-id", po::value<std::string>(&run_id)->default_value("run"),
      "Structured output run ID")(
      "output-dir",
      po::value<std::string>(&output_dir)->default_value("results"),
      "Structured output root directory")(
      "passes,p", po::value<long long>(&passes)->default_value(100),
      "Total pass target (resume default: saved target)")(
      "checkpoint,c", po::value<long long>(&checkpoint)->default_value(10),
      "Checkpoint every n global passes");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, description), args);

  if (args.count("help"))
  {
    fmt::print("{}\n", fmt::streamed(description));
    return EXIT_SUCCESS;
  }

  if (args.count("version"))
  {
    fmt::print("CDT++ version {}\n", cdt::VERSION);
    return EXIT_SUCCESS;
  }

  po::notify(args);
  auto const has_input           = args.count("input") != 0;
  auto const has_resume          = args.count("resume") != 0;
  auto const explicitly_supplied = [&args](char const* option) {
    auto const value = args.find(option);
    return value != args.end() && !value->second.defaulted();
  };
  if (has_input && has_resume)
  {
    throw invalid_argument("--input and --resume are mutually exclusive.");
  }
  if ((has_input || has_resume) &&
      (args.count("spherical") != 0 || args.count("toroidal") != 0 ||
       args.count("simplices") != 0 || args.count("timeslices") != 0 ||
       explicitly_supplied("dimensions") || explicitly_supplied("init") ||
       explicitly_supplied("foliate")))
  {
    throw invalid_argument(fmt::format(
        "{} cannot be combined with topology or triangulation-construction "
        "options.",
        has_resume ? "--resume" : "--input"));
  }
  if (has_resume &&
      (args.count("seed") != 0 || explicitly_supplied("threads") ||
       args.count("alpha") != 0 || args.count("k") != 0 ||
       args.count("lambda") != 0 || args.count("kappa0") != 0 ||
       args.count("kappa4") != 0 || args.count("Delta") != 0 ||
       explicitly_supplied("target-n4") ||
       explicitly_supplied("volume-epsilon") ||
       explicitly_supplied("thermalization") ||
       explicitly_supplied("measurement-interval") ||
       explicitly_supplied("checkpoint")))
  {
    throw invalid_argument(
        "--resume restores seed, threads, action parameters, and checkpoint "
        "cadence from the saved run.");
  }
  if (dimensions != 3 && dimensions != 4)
  {
    throw invalid_argument(
        "Only three- or four-dimensional triangulations are supported.");
  }
  if (dimensions == 4)
  {
    if (has_input || has_resume)
    {
      throw invalid_argument("4D runs do not support --input or --resume yet.");
    }
    if (!args.count("simplices"))
    {
      throw invalid_argument("Number of simplices not specified.");
    }
    if (!args.count("timeslices"))
    {
      throw invalid_argument("Number of timeslices not specified.");
    }
    if (args.count("kappa0") == 0 || args.count("kappa4") == 0 ||
        args.count("Delta") == 0)
    {
      throw invalid_argument(
          "4D runs require explicit --kappa0, --kappa4, and --Delta.");
    }
    auto root_random =
        args.count("seed") != 0 ? cdt::Random{seed} : cdt::Random{};
    auto const write_files = !args.count("no-output");
    auto const config = runtime_config::make_4d_simulation(
        args.count("spherical") != 0, args.count("toroidal") != 0,
        simplices, timeslices, dimensions, kappa_0, kappa_4, Delta,
        target_N4, volume_epsilon, passes, checkpoint, thermalization,
        measurement_interval, root_random.seed(), chain_id, run_id,
        std::filesystem::path{output_dir}, threads, write_files);

    fmt::print("Topology is spherical\n");
    fmt::print("Dimensionality: 3+1\n");
    fmt::print("Number of desired simplices: {}\n", config.simplices());
    fmt::print("Number of desired timeslices: {}\n", config.timeslices());
    fmt::print("Number of passes: {}\n", config.passes());
    fmt::print("4D sweep size: {} move attempts.\n", config.sweep_size());
    fmt::print("Total 4D move attempts: {}\n", config.steps());
    fmt::print("Checkpoint every {} move attempts.\n",
               config.checkpoint_steps());
    fmt::print("Effective random seed: {}\n", root_random.seed());
    fmt::print("=== Parameters ===\n");
    fmt::print("kappa_0: {}\n", config.metropolis().couplings.kappa_0);
    fmt::print("kappa_4: {}\n", config.metropolis().couplings.kappa_4);
    fmt::print("Delta: {}\n", config.metropolis().couplings.Delta);

    Timer timer;
    timer.start();
    fmt::print("cdt started at {}\n", utilities::current_date_time());

    auto universe =
        four_d::FoliatedTriangulation4::periodic_seed(config.timeslices());

    four_d::Metropolis4 run(config.metropolis());
    auto result = run.run(std::move(universe), config.steps());
    if (!result.triangulation.is_valid())
    {
      throw runtime_error("4D result is invalid!\n");
    }

    if (config.write_files())
    {
      auto const run_dir = config.output_dir() / config.run_id();
      run.save_checkpoint(run_dir / "checkpoint", result.triangulation,
                          config.steps());
      four_d::output::RunManifest manifest;
      manifest.run_id     = config.run_id();
      manifest.git_commit = std::string(cdt::SOURCE_REVISION);
      manifest.build_type = std::string(cdt::BUILD_CONFIGURATION);
      manifest.compiler =
          fmt::format("{} {}", cdt::BUILD_COMPILER_ID,
                      cdt::BUILD_COMPILER_VERSION);
      four_d::output::write_run_directory(config.output_dir(), manifest,
                                          config.metropolis(), result);
    }

    timer.stop();
    fmt::print("=== 4D Run Results ===\n");
    fmt::print("Running time is {} seconds.\n", timer.time());
    if (config.write_files())
    {
      fmt::print("Structured output written to {}/{}\n",
                 config.output_dir().string(), config.run_id());
    }
    else
    {
      fmt::print("Structured output disabled.\n");
    }
    Int_precision attempted{0};
    Int_precision accepted{0};
    Int_precision invalid{0};
    for (auto const& stat : result.move_stats)
    {
      attempted += stat.attempted;
      accepted += stat.accepted;
      invalid += stat.invalid;
    }
    fmt::print("Move attempts: {}\n", attempted);
    fmt::print("Accepted moves: {}\n", accepted);
    fmt::print("Invalid proposals: {}\n", invalid);
    fmt::print("Final N4: {}\n", result.triangulation.counts().N4);
    fmt::print("Measurements: {}\n", result.measurements.size());
    auto const report = result.triangulation.validate();
    fmt::print("Standard CDT candidate: {}\n",
               report.valid() && report.standard_cdt_candidate ? "true"
                                                               : "false");
    return EXIT_SUCCESS;
  }
  if (!has_resume && (args.count("alpha") == 0 || args.count("k") == 0 ||
                      args.count("lambda") == 0))
  {
    throw invalid_argument("Alpha, K, and Lambda must be specified.");
  }
  if (!has_input && !has_resume && !args.count("simplices"))
  {
    throw invalid_argument("Number of simplices not specified.");
  }
  if (!has_input && !has_resume && !args.count("timeslices"))
  {
    throw invalid_argument("Number of timeslices not specified.");
  }

  using Initial_artifact =
      utilities::Initial_triangulation_artifact<Delaunay_t<3>>;
  using Resume_artifact = utilities::Checkpoint_artifact<Delaunay_t<3>>;
  std::optional<Initial_artifact> initial_artifact;
  std::optional<Resume_artifact>  resume_artifact;
  if (has_input)
  {
    initial_artifact.emplace(
        utilities::read_initial_triangulation<Delaunay_t<3>>(input_path));
  }
  if (has_resume)
  {
    resume_artifact.emplace(
        utilities::read_checkpoint<Delaunay_t<3>>(resume_path));
  }

  auto root_random =
      resume_artifact
          ? cdt::Random{resume_artifact->metadata.seed}
          : (args.count("seed") != 0 ? cdt::Random{seed} : cdt::Random{});
  auto const effective_threads = [&] {
    if (!resume_artifact) { return threads; }
    auto const saved = *resume_artifact->metadata.max_threads;
    if (!std::in_range<long long>(saved))
    {
      throw out_of_range("Saved thread count exceeds the supported range.");
    }
    return static_cast<long long>(saved);
  }();
  auto const triangulation_config = [&] {
    if (initial_artifact || resume_artifact)
    {
      auto const& metadata = initial_artifact ? initial_artifact->metadata
                                              : resume_artifact->metadata;
      return runtime_config::make_triangulation(
          metadata.topology == Topology::SPHERICAL,
          metadata.topology == Topology::TOROIDAL, metadata.desired_simplices,
          metadata.desired_timeslices, metadata.dimension,
          metadata.initial_radius, metadata.foliation_spacing,
          root_random.seed(), effective_threads);
    }
    return runtime_config::make_triangulation(
        args.count("spherical") != 0, args.count("toroidal") != 0, simplices,
        timeslices, dimensions, initial_radius, foliation_spacing,
        root_random.seed(), effective_threads);
  }();
  auto const effective_alpha =
      resume_artifact ? *resume_artifact->metadata.alpha : alpha;
  auto const effective_k = resume_artifact ? *resume_artifact->metadata.k : k;
  auto const effective_lambda =
      resume_artifact ? *resume_artifact->metadata.lambda : lambda;
  auto const completed_passes = resume_artifact
                                  ? *resume_artifact->metadata.completed_passes
                                  : Int_precision{};
  auto const target_passes    = [&] {
    if (!resume_artifact) { return passes; }
    if (explicitly_supplied("passes")) { return passes; }
    return static_cast<long long>(*resume_artifact->metadata.configured_passes);
  }();
  if (resume_artifact && (!std::in_range<Int_precision>(target_passes) ||
                          target_passes < completed_passes))
  {
    throw invalid_argument(
        "Resume target passes must be at least the completed checkpoint pass.");
  }
  auto const passes_to_execute =
      target_passes - static_cast<long long>(completed_passes);
  auto const effective_checkpoint =
      resume_artifact ? static_cast<long long>(
                            *resume_artifact->metadata.checkpoint_interval)
                      : checkpoint;
  auto const config = runtime_config::make_simulation(
      triangulation_config, effective_alpha, effective_k, effective_lambda,
      resume_artifact && passes_to_execute == 0 ? 1 : passes_to_execute,
      effective_checkpoint, !args.count("no-output"));
#if defined(CDT_ENABLE_PARALLEL_TRIANGULATION) && \
    CDT_ENABLE_PARALLEL_TRIANGULATION
  [[maybe_unused]] oneapi::tbb::global_control thread_limit{
      oneapi::tbb::global_control::max_allowed_parallelism,
      config.triangulation().threads()};
#endif
  auto transition_random =
      resume_artifact ? cdt::Random::from_serialized_state(
                            resume_artifact->metadata.seed,
                            resume_artifact->metadata.transition_stream,
                            *resume_artifact->metadata.transition_random_state)
                      : root_random.split(cdt::random_streams::transitions);

  // Display job parameters
  fmt::print("Topology is {}\n",
             utilities::topology_to_str(config.triangulation().topology()));
  fmt::print("Dimensionality: {}+{}\n", config.triangulation().dimensions() - 1,
             1);
  fmt::print("Initial radius: {}\n", config.triangulation().initial_radius());
  fmt::print("Foliation spacing: {}\n",
             config.triangulation().foliation_spacing());
  fmt::print("Number of desired simplices: {}\n",
             config.triangulation().simplices());
  fmt::print("Number of desired timeslices: {}\n",
             config.triangulation().timeslices());
  fmt::print("Number of passes to execute: {}\n", passes_to_execute);
  fmt::print("Checkpoint every {} passes.\n", config.checkpoint());
  fmt::print("Effective random seed: {}\n", config.triangulation().seed());
  fmt::print("Maximum Delaunay threads: {}\n",
             config.triangulation().threads());
  if (initial_artifact)
  {
    fmt::print("Input initial triangulation: {}\n", input_path);
    fmt::print("Input initialization seed: {}\n",
               initial_artifact->metadata.seed);
    fmt::print("Input topology fingerprint: {:016x}\n",
               *initial_artifact->metadata.topology_fingerprint);
  }
  if (resume_artifact)
  {
    fmt::print("Resuming checkpoint: {}\n", resume_path);
    fmt::print("Completed checkpoint passes: {}\n", completed_passes);
    fmt::print("Target total passes: {}\n", target_passes);
    fmt::print("Checkpoint transition count: {}\n",
               *resume_artifact->metadata.transition_count);
  }
  fmt::print("=== Parameters ===\n");
  fmt::print("Alpha: {}\n", config.alpha());
  fmt::print("K: {}\n", config.k());
  fmt::print("Lambda: {}\n", config.lambda());

  // Start running time
  Timer timer;
  timer.start();
  fmt::print("cdt started at {}\n", utilities::current_date_time());

  // Load an exact checkpoint, load an initial state, or generate a fresh state.
  auto universe = [&]() -> manifolds::Manifold_3 {
    if (initial_artifact || resume_artifact)
    {
      auto const& metadata      = initial_artifact ? initial_artifact->metadata
                                                   : resume_artifact->metadata;
      auto        triangulation = initial_artifact
                                    ? std::move(initial_artifact->triangulation)
                                    : std::move(resume_artifact->triangulation);
      auto        manifold      = manifolds::Manifold_3{
          foliated_triangulations::FoliatedTriangulation_3{
                                                           std::move(triangulation), metadata.initial_radius,
                                                           metadata.foliation_spacing}
      };
      if (!manifold.is_correct_with_diagnostics())
      {
        throw invalid_argument(
            "Input triangulation does not satisfy the CDT manifold contract.");
      }
      return manifold;
    }
    auto initialization_random =
        root_random.split(cdt::random_streams::initialization);
    return manifolds::Manifold_3{
        config.triangulation().simplices(), config.triangulation().timeslices(),
        initialization_random, config.triangulation().initial_radius(),
        config.triangulation().foliation_spacing()};
  }();

  auto reproducibility     = resume_artifact
                               ? resume_artifact->metadata
                               : utilities::make_reproducibility_metadata(
                                     universe, config.triangulation().seed(),
                                     utilities::ArtifactKind::FINAL_TRIANGULATION);
  reproducibility.artifact = utilities::ArtifactKind::FINAL_TRIANGULATION;
  reproducibility.desired_simplices  = config.triangulation().simplices();
  reproducibility.desired_timeslices = config.triangulation().timeslices();
  reproducibility.alpha              = config.alpha();
  reproducibility.k                  = config.k();
  reproducibility.lambda             = config.lambda();
  reproducibility.configured_passes = static_cast<Int_precision>(target_passes);
  reproducibility.checkpoint_interval = config.checkpoint();
  reproducibility.max_threads         = config.triangulation().threads();
  if (!resume_artifact)
  {
    reproducibility.input_artifact =
        utilities::ArtifactKind::INITIAL_TRIANGULATION;
    reproducibility.input_seed = initial_artifact
                                   ? initial_artifact->metadata.seed
                                   : config.triangulation().seed();
    reproducibility.input_initialization_stream =
        initial_artifact ? initial_artifact->metadata.initialization_stream
                         : cdt::random_streams::initialization;
    reproducibility.input_placement_fingerprint =
        initial_artifact ? initial_artifact->metadata.placement_fingerprint
                         : reproducibility.placement_fingerprint;
    reproducibility.input_topology_fingerprint =
        initial_artifact ? initial_artifact->metadata.topology_fingerprint
                         : reproducibility.topology_fingerprint;
  }

  // Look at triangulation
  universe.print();
  universe.print_details();
  universe.print_volume_per_timeslice();

  if (resume_artifact && passes_to_execute == 0)
  {
    fmt::print(
        "Checkpoint already reached the target pass; no transitions remain.\n");
    if (config.write_files())
    {
      reproducibility.transition_random_state.reset();
      utilities::write_file(universe, reproducibility);
    }
    return EXIT_SUCCESS;
  }

  // Initialize the Metropolis algorithm with complete run provenance.
  Metropolis_3 run(config.alpha(), config.k(), config.lambda(), config.passes(),
                   config.checkpoint(), config.write_files(),
                   std::move(transition_random), reproducibility,
                   completed_passes);

  // The main work of the program
  auto const result = run(universe);

  // Do we have enough timeslices?
  if (auto max_timevalue = result.max_time();
      max_timevalue < config.triangulation().timeslices())
  {
    fmt::print("You wanted {} timeslices, but only got {}.\n",
               config.triangulation().timeslices(), max_timevalue);
  }

  if (!result.is_valid()) { throw runtime_error("Result is invalid!\n"); }

  // Print results
  timer.stop();  // End running time counter
  fmt::print("=== Run Results ===\n");
  fmt::print("Running time is {} seconds.\n", timer.time());
  result.print();
  result.print_details();
  result.print_volume_per_timeslice();

  // Write results to file
  if (config.write_files())
  {
    utilities::write_file(
        result, run.reproducibility_metadata(
                    result, utilities::ArtifactKind::FINAL_TRIANGULATION,
                    static_cast<Int_precision>(target_passes)));
  }

  return EXIT_SUCCESS;
}
catch (domain_error const& DomainError)
{
  spdlog::critical("{}\n", DomainError.what());
  spdlog::critical("Triangle inequalities violated ... Exiting.\n");
  return EXIT_FAILURE;
}
catch (invalid_argument const& InvalidArgument)
{
  spdlog::critical("{}\n", InvalidArgument.what());
  spdlog::critical("Invalid parameter ... Exiting.\n");
  return EXIT_FAILURE;
}
catch (logic_error const& LogicError)
{
  spdlog::critical("{}\n", LogicError.what());
  spdlog::critical("Simulation startup failed ... Exiting.\n");
  return EXIT_FAILURE;
}
catch (runtime_error const& RuntimeError)
{
  spdlog::critical("{}\n", RuntimeError.what());
  return EXIT_FAILURE;
}
catch (...)
{
  spdlog::critical("Something went wrong ... Exiting.\n");
  return EXIT_FAILURE;
}
