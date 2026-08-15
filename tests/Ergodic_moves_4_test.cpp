#include "Ergodic_moves_4.hpp"

#include <doctest/doctest.h>

using namespace cdt::four_d;
namespace move_tracker = cdt::move_tracker;

namespace
{
  void check_reverse_site_restores(FoliatedTriangulation4 const& before,
                                   moves::MoveApplication const& moved)
  {
    auto const reverse = move_tracker::reverse_move(moved.move);
    auto const reverse_candidates =
        moved.triangulation.candidate_multiplicity(reverse);
    auto restored = false;
    for (auto site = std::size_t{0};
         site < static_cast<std::size_t>(reverse_candidates); ++site)
    {
      auto inverse = moves::apply(moved.triangulation, reverse, site);
      if (inverse &&
          inverse->triangulation.canonical_hash() == before.canonical_hash())
      {
        restored = inverse->triangulation.is_valid();
        break;
      }
    }
    CHECK(restored);
  }

  void check_delta_matches_counts(FoliatedTriangulation4 const& before,
                                  moves::MoveApplication const& moved)
  {
    auto const before_counts = before.counts();
    auto const after_counts  = moved.triangulation.counts();
    CHECK_EQ(moved.delta.N0, after_counts.N0 - before_counts.N0);
    CHECK_EQ(moved.delta.N1, after_counts.N1 - before_counts.N1);
    CHECK_EQ(moved.delta.N2, after_counts.N2 - before_counts.N2);
    CHECK_EQ(moved.delta.N3, after_counts.N3 - before_counts.N3);
    CHECK_EQ(moved.delta.N4, after_counts.N4 - before_counts.N4);
    CHECK_EQ(moved.delta.N41, after_counts.N41 - before_counts.N41);
    CHECK_EQ(moved.delta.N32, after_counts.N32 - before_counts.N32);
    CHECK_EQ(moved.delta.N23, after_counts.N23 - before_counts.N23);
    CHECK_EQ(moved.delta.N14, after_counts.N14 - before_counts.N14);
  }
}  // namespace

TEST_CASE("4D reverse_move covers every inverse pair")
{
  using move_tracker::MoveType4D;
  CHECK(move_tracker::reverse_move(MoveType4D::TWO_FOUR) ==
        MoveType4D::FOUR_TWO);
  CHECK(move_tracker::reverse_move(MoveType4D::FOUR_TWO) ==
        MoveType4D::TWO_FOUR);
  CHECK(move_tracker::reverse_move(MoveType4D::THREE_THREE) ==
        MoveType4D::THREE_THREE);
  CHECK(move_tracker::reverse_move(MoveType4D::FOUR_SIX) ==
        MoveType4D::SIX_FOUR);
  CHECK(move_tracker::reverse_move(MoveType4D::SIX_FOUR) ==
        MoveType4D::FOUR_SIX);
  CHECK(move_tracker::reverse_move(MoveType4D::TWO_EIGHT) ==
        MoveType4D::EIGHT_TWO);
  CHECK(move_tracker::reverse_move(MoveType4D::EIGHT_TWO) ==
        MoveType4D::TWO_EIGHT);
}

TEST_CASE("4D move plus inverse restores canonical hash")
{
  using move_tracker::MoveType4D;
  auto seed     = FoliatedTriangulation4::periodic_seed(4);
  auto two_four = moves::apply(seed, MoveType4D::TWO_FOUR);
  REQUIRE(two_four);
  CHECK(two_four->triangulation.is_valid());
  check_reverse_site_restores(seed, *two_four);

  auto four_two = moves::apply(two_four->triangulation, MoveType4D::FOUR_TWO);
  REQUIRE(four_two);
  CHECK(four_two->triangulation.is_valid());
  check_reverse_site_restores(two_four->triangulation, *four_two);
}

TEST_CASE("4D moves have exact incidence-derived combinatorial count changes")
{
  auto seed     = FoliatedTriangulation4::periodic_seed(4);
  auto two_four = moves::apply(seed, move_tracker::MoveType4D::TWO_FOUR);
  REQUIRE(two_four);
  check_delta_matches_counts(seed, *two_four);
  CHECK_EQ(two_four->delta.N0, 0);
  CHECK_EQ(two_four->delta.N1, 1);
  CHECK_EQ(two_four->delta.N2, 4);
  CHECK_EQ(two_four->delta.N3, 5);
  CHECK_EQ(two_four->delta.N4, 2);
  CHECK_EQ(two_four->delta.N41 + two_four->delta.N32 + two_four->delta.N23 +
               two_four->delta.N14,
           two_four->delta.N4);

  auto four_two =
      moves::apply(two_four->triangulation, move_tracker::MoveType4D::FOUR_TWO);
  REQUIRE(four_two);
  check_delta_matches_counts(two_four->triangulation, *four_two);
  CHECK_EQ(four_two->delta.N0, 0);
  CHECK_EQ(four_two->delta.N1, -1);
  CHECK_EQ(four_two->delta.N2, -4);
  CHECK_EQ(four_two->delta.N3, -5);
  CHECK_EQ(four_two->delta.N4, -2);
  CHECK_EQ(four_two->delta.N41 + four_two->delta.N32 + four_two->delta.N23 +
               four_two->delta.N14,
           four_two->delta.N4);
}

TEST_CASE("4D local and full action differences agree for every move")
{
  auto        seed = FoliatedTriangulation4::periodic_seed(4);
  S4Couplings couplings{1.0L, 0.2L, 0.1L, 64, 0.001L};
  auto        check_action_delta = [&](FoliatedTriangulation4 const& before,
                                moves::MoveApplication const& moved) {
    auto const full_delta = S4_action_difference(
        before.counts(), moved.triangulation.counts(), couplings);
    auto const local_delta =
        local_action_difference(before.counts(), moved.delta, couplings);
    CHECK(full_delta == doctest::Approx(local_delta));
  };

  auto two_four = moves::apply(seed, move_tracker::MoveType4D::TWO_FOUR);
  REQUIRE(two_four);
  check_action_delta(seed, *two_four);

  auto four_two =
      moves::apply(two_four->triangulation, move_tracker::MoveType4D::FOUR_TWO);
  REQUIRE(four_two);
  check_action_delta(two_four->triangulation, *four_two);
}

TEST_CASE("4D forward and reverse proposal multiplicities are state-derived")
{
  auto       seed = FoliatedTriangulation4::periodic_seed(4);
  auto const forward =
      seed.candidate_multiplicity(move_tracker::MoveType4D::TWO_FOUR);
  auto two_four = moves::apply(seed, move_tracker::MoveType4D::TWO_FOUR);
  REQUIRE(two_four);
  auto const reverse = two_four->triangulation.candidate_multiplicity(
      move_tracker::MoveType4D::FOUR_TWO);
  CHECK_EQ(two_four->forward_candidates, forward);
  CHECK_EQ(two_four->reverse_candidates, reverse);
  CHECK_GT(forward, 0);
  CHECK_EQ(reverse, 1);
}

TEST_CASE("Unsupported 4D move descriptors have no production proposal sites")
{
  auto triangulation = FoliatedTriangulation4::periodic_seed(4);
  using move_tracker::MoveType4D;
  for (auto const move :
       {MoveType4D::THREE_THREE, MoveType4D::FOUR_SIX, MoveType4D::SIX_FOUR,
        MoveType4D::TWO_EIGHT, MoveType4D::EIGHT_TWO})
  {
    CHECK_EQ(triangulation.candidate_multiplicity(move), 0);
    CHECK_FALSE(moves::apply(triangulation, move));
  }
}

TEST_CASE("4D failed moves leave the original unchanged")
{
  auto                   counts = S4Counts{1, 0, 0, 0, 0, 0, 0, 0, 0};
  FoliatedTriangulation4 triangulation{
      4, counts, {1, 1, 1, 1}
  };
  auto const hash = triangulation.canonical_hash();
  auto result = moves::apply(triangulation, move_tracker::MoveType4D::FOUR_TWO);
  CHECK_FALSE(result);
  CHECK_EQ(triangulation.canonical_hash(), hash);
}
