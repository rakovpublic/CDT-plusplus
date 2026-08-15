#include "Foliated_triangulation_4.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <string>
#include <string_view>

using namespace cdt::four_d;
namespace move_tracker = cdt::move_tracker;

namespace
{
  [[nodiscard]] auto has_error(ValidationReport const& report,
                               std::string_view const  error) -> bool
  {
    return std::ranges::any_of(
        report.errors,
        [error](std::string const& candidate) { return candidate == error; });
  }

  [[nodiscard]] auto isolated_simplex(
      SimplexId const id, std::array<VertexId, 5> vertices,
      SimplexType4D const type = SimplexType4D::FOUR_ONE) -> Simplex4D
  {
    Simplex4D simplex;
    simplex.id       = id;
    simplex.vertices = vertices;
    simplex.type     = type;
    return simplex;
  }
}  // namespace

TEST_CASE("Persistent 4D periodic S3xS1 seed validates")
{
  auto triangulation = FoliatedTriangulation4::periodic_seed(4);
  auto counts        = triangulation.counts();

  CHECK(triangulation.periodic());
  CHECK_EQ(triangulation.timeslices(), 4);
  CHECK_FALSE(triangulation.vertices().empty());
  CHECK_FALSE(triangulation.simplices().empty());
  CHECK_EQ(triangulation.spatial_topology(), "S3");
  CHECK_EQ(triangulation.spacetime_topology(), "S3xS1");
  CHECK(triangulation.has_closed_s3_slices());
  CHECK(triangulation.is_valid());
  CHECK_EQ(counts.N0, 20);
  CHECK_EQ(counts.N4, 80);
  CHECK_EQ(counts.N4, counts.N41 + counts.N32 + counts.N23 + counts.N14);
  CHECK_GT(counts.N41, 0);
  CHECK_GT(counts.N32, 0);
  CHECK_GT(counts.N23, 0);
  CHECK_GT(counts.N14, 0);
  CHECK_EQ(triangulation.occupied_temporal_width(), 4);
  for (auto const chi : triangulation.slice_euler_characteristics())
  {
    CHECK_EQ(chi, 0);
  }
}

TEST_CASE("Persistent 4D canonical hash tracks local incidence changes")
{
  auto triangulation = FoliatedTriangulation4::periodic_seed(3);
  auto copy          = triangulation;
  CHECK_EQ(copy.canonical_hash(), triangulation.canonical_hash());

  auto moved = triangulation;
  REQUIRE(moved.apply_move(move_tracker::MoveType4D::TWO_FOUR));
  CHECK_NE(moved.canonical_hash(), triangulation.canonical_hash());
  auto const before_counts = triangulation.counts();
  auto const after_counts  = moved.counts();
  CHECK_EQ(after_counts.N0, before_counts.N0);
  CHECK_EQ(after_counts.N1, before_counts.N1 + 1);
  CHECK_EQ(after_counts.N4, before_counts.N4 + 2);
  CHECK_FALSE(moved.simplices().empty());
}

TEST_CASE("Unsupported 4D moves are not advertised as local proposals")
{
  auto triangulation = FoliatedTriangulation4::periodic_seed(3);
  CHECK_EQ(triangulation.candidate_multiplicity(
               move_tracker::MoveType4D::THREE_THREE),
           0);
  CHECK_EQ(
      triangulation.candidate_multiplicity(move_tracker::MoveType4D::TWO_EIGHT),
      0);
  CHECK_FALSE(triangulation.apply_move(move_tracker::MoveType4D::THREE_THREE));
  CHECK_FALSE(triangulation.apply_move(move_tracker::MoveType4D::TWO_EIGHT));
}

TEST_CASE("4D time reversal maps profiles and vertex times cyclically")
{
  auto       triangulation = FoliatedTriangulation4::periodic_seed(4);
  auto       reversed      = triangulation.time_reversed();
  auto const profile       = reversed.spatial_volume_profile();
  auto const seed_profile  = triangulation.spatial_volume_profile();
  REQUIRE_EQ(profile.size(), 4);
  CHECK_EQ(profile[0], seed_profile[0]);
  CHECK_EQ(profile[1], seed_profile[3]);
  CHECK_EQ(profile[2], seed_profile[2]);
  CHECK_EQ(profile[3], seed_profile[1]);
  REQUIRE_EQ(reversed.vertices().size(), triangulation.vertices().size());
  CHECK_EQ(reversed.vertices()[0].time, 0);
  CHECK_EQ(reversed.vertices()[5].time, 3);
  CHECK_EQ(reversed.vertices()[10].time, 2);
  CHECK_EQ(reversed.vertices()[15].time, 1);
  CHECK_FALSE(reversed.three_three_forward());
  CHECK(reversed.is_valid());
}

TEST_CASE("4D candidate validation is independent from the initializer")
{
  auto seeded = FoliatedTriangulation4::periodic_seed(3);

  SUBCASE("count-only states are not standard CDT candidates")
  {
    auto from_counts = FoliatedTriangulation4::from_counts_for_validation(
        seeded.timeslices(), seeded.counts(), seeded.spatial_volume_profile());

    auto const report = from_counts.validate();
    CHECK_FALSE(report.valid());
    CHECK_FALSE(report.standard_cdt_candidate);
    CHECK_EQ(from_counts.spatial_topology(), "non-S3");
    CHECK_EQ(from_counts.spacetime_topology(), "unvalidated");
    CHECK_EQ(
        from_counts.candidate_multiplicity(move_tracker::MoveType4D::TWO_FOUR),
        0);
  }

  SUBCASE("count-only constructor clamps timeslices")
  {
    auto clamped = FoliatedTriangulation4::from_counts_for_validation(
        -7, seeded.counts(), {});
    CHECK_EQ(clamped.timeslices(), 2);
    CHECK_EQ(clamped.spatial_volume_profile().size(), 2);
  }

  SUBCASE("count-only proposal inventories do not infer local move sites")
  {
    auto counts   = S4Counts{1, 2, 3, 4, 5, 1, 2, 1, 1};
    auto abstract = FoliatedTriangulation4::from_counts_for_validation(
        2, counts, FoliatedTriangulation4::Profile{1, 1});
    auto const inventory = abstract.proposal_inventory();
    CHECK_EQ(inventory.spatial_tetrahedra, 0);
    CHECK_EQ(inventory.timelike_edges, 0);
    CHECK_EQ(inventory.mixed_triangles, 0);
    CHECK_EQ(inventory.timelike_tetrahedra, 0);
    CHECK_EQ(inventory.vertices, 0);
    CHECK_EQ(inventory.three_two_simplices, 0);
    CHECK_EQ(inventory.two_three_simplices, 0);
  }

  SUBCASE("class-resolved counts still do not replace local site enumeration")
  {
    auto counts           = S4Counts{1, 2, 3, 4, 5, 1, 2, 1, 1};
    counts.class_resolved = S4ClassResolvedCounts{7, 8, 9, 10};
    auto exact            = FoliatedTriangulation4::from_counts_for_validation(
        2, counts, FoliatedTriangulation4::Profile{1, 1});
    auto const exact_inventory = exact.proposal_inventory();
    CHECK_EQ(exact_inventory.spatial_tetrahedra, 0);
    CHECK_EQ(exact_inventory.timelike_edges, 0);
    CHECK_EQ(exact_inventory.mixed_triangles, 0);
    CHECK_EQ(exact_inventory.timelike_tetrahedra, 0);
  }

  SUBCASE("negative spatial profile is reported")
  {
    auto invalid = FoliatedTriangulation4::from_counts_for_validation(
        2, seeded.counts(), FoliatedTriangulation4::Profile{-1, 1});
    auto const report = invalid.validate();
    CHECK_FALSE(report.valid());
    CHECK(has_error(report, "Spatial profile contains negative volume."));
  }

  SUBCASE("vertex-time cache errors are reported")
  {
    auto const empty_counts = S4Counts{};
    auto       invalid      = FoliatedTriangulation4::from_checkpoint_state(
        2, empty_counts,
        FoliatedTriangulation4::Profile{
            0, 0
    },
        FoliatedTriangulation4::VertexContainer{Vertex4D{1, 0}, Vertex4D{1, 1}},
        {}, true);
    auto const report = invalid.validate();
    CHECK_FALSE(report.valid());
    CHECK(has_error(report, "Vertex-time cache does not match vertices."));
    CHECK(has_error(report, "Vertex-time cache is stale."));
  }

  SUBCASE("missing simplex vertices are reported")
  {
    auto const empty_counts = S4Counts{};
    auto       invalid      = FoliatedTriangulation4::from_checkpoint_state(
        2, empty_counts,
        FoliatedTriangulation4::Profile{
            0, 0
    },
        FoliatedTriangulation4::VertexContainer{Vertex4D{1, 0}},
        FoliatedTriangulation4::SimplexContainer{
            isolated_simplex(1, std::array<VertexId, 5>{1, 2, 3, 4, 5})},
        true);

    auto const report = invalid.validate();

    CHECK_FALSE(report.valid());
    CHECK(has_error(report, "A 4-simplex references a missing vertex."));
  }

  SUBCASE("disconnected restored complexes are rejected")
  {
    auto const counts   = S4Counts{999, 0, 0, 0, 999, 0, 0, 0, 0};
    auto const vertices = FoliatedTriangulation4::VertexContainer{
        Vertex4D{ 1, 0},
        Vertex4D{ 2, 0},
        Vertex4D{ 3, 0},
        Vertex4D{ 4, 0},
        Vertex4D{ 5, 1},
        Vertex4D{ 6, 0},
        Vertex4D{ 7, 0},
        Vertex4D{ 8, 0},
        Vertex4D{ 9, 0},
        Vertex4D{10, 1}
    };
    auto invalid = FoliatedTriangulation4::from_checkpoint_state(
        2, counts, FoliatedTriangulation4::Profile{2, 0}, vertices,
        FoliatedTriangulation4::SimplexContainer{
            isolated_simplex(1, std::array<VertexId, 5>{1, 2, 3, 4, 5}),
            isolated_simplex(2, std::array<VertexId, 5>{6, 7, 8, 9, 10})},
        true);
    auto const invalid_counts = invalid.counts();
    CHECK_EQ(invalid_counts.N0, 10);
    CHECK_EQ(invalid_counts.N4, 2);
    CHECK_EQ(invalid_counts.N41, 2);
    REQUIRE(invalid_counts.class_resolved.has_value());
    CHECK_EQ(invalid.proposal_inventory().spatial_tetrahedra, 0);

    auto const report = invalid.validate();
    CHECK_FALSE(report.valid());
    CHECK(has_error(report, "Simplex neighbor graph is disconnected."));
    CHECK(has_error(
        report, "Spatial slices are not validated as connected S3 slices."));
  }
}
