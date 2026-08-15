#include "Ergodic_moves_4.hpp"

#include <doctest/doctest.h>

#include <optional>

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

  [[nodiscard]] auto state_with_four_six_site(
      FoliatedTriangulation4 const& seed)
      -> std::optional<FoliatedTriangulation4>
  {
    auto const first_sites =
        seed.candidate_multiplicity(move_tracker::MoveType4D::TWO_EIGHT);
    for (auto first_site = std::size_t{0};
         first_site < static_cast<std::size_t>(first_sites); ++first_site)
    {
      auto first = moves::apply(seed, move_tracker::MoveType4D::TWO_EIGHT,
                                first_site);
      if (!first) { continue; }
      auto const second_sites = first->triangulation.candidate_multiplicity(
          move_tracker::MoveType4D::TWO_EIGHT);
      for (auto second_site = std::size_t{0};
           second_site < static_cast<std::size_t>(second_sites); ++second_site)
      {
        auto second =
            moves::apply(first->triangulation,
                         move_tracker::MoveType4D::TWO_EIGHT, second_site);
        if (second && second->triangulation.candidate_multiplicity(
                          move_tracker::MoveType4D::FOUR_SIX) > 0)
        {
          return second->triangulation;
        }
      }
    }
    return std::nullopt;
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
  auto seed = FoliatedTriangulation4::periodic_seed(4);
  auto two_four = moves::apply(seed, MoveType4D::TWO_FOUR);
  REQUIRE(two_four);
  CHECK(two_four->triangulation.is_valid());
  check_reverse_site_restores(seed, *two_four);

  auto three_three = moves::apply(two_four->triangulation,
                                  MoveType4D::THREE_THREE);
  REQUIRE(three_three);
  CHECK(three_three->triangulation.is_valid());
  check_reverse_site_restores(two_four->triangulation, *three_three);

  auto two_eight = moves::apply(seed, MoveType4D::TWO_EIGHT);
  REQUIRE(two_eight);
  CHECK(two_eight->triangulation.is_valid());
  check_reverse_site_restores(seed, *two_eight);

  auto four_six_before = state_with_four_six_site(seed);
  REQUIRE(four_six_before);
  auto four_six = moves::apply(*four_six_before, MoveType4D::FOUR_SIX);
  REQUIRE(four_six);
  CHECK(four_six->triangulation.is_valid());
  check_reverse_site_restores(*four_six_before, *four_six);
}

TEST_CASE("4D moves have exact incidence-derived combinatorial count changes")
{
  auto seed = FoliatedTriangulation4::periodic_seed(4);
  auto check_expected_delta = [&](FoliatedTriangulation4 const& before,
                                  move_tracker::MoveType4D const move,
                                  S4Counts const& expected_delta) {
    auto moved = moves::apply(before, move);
    REQUIRE(moved);
    check_delta_matches_counts(before, *moved);
    CHECK_EQ(moved->delta.N0, expected_delta.N0);
    CHECK_EQ(moved->delta.N1, expected_delta.N1);
    CHECK_EQ(moved->delta.N2, expected_delta.N2);
    CHECK_EQ(moved->delta.N3, expected_delta.N3);
    CHECK_EQ(moved->delta.N4, expected_delta.N4);
    CHECK_EQ(moved->delta.N41 + moved->delta.N32 + moved->delta.N23 +
                 moved->delta.N14,
             moved->delta.N4);

    auto inverse = moves::apply(moved->triangulation,
                                move_tracker::reverse_move(move));
    REQUIRE(inverse);
    check_delta_matches_counts(moved->triangulation, *inverse);
    CHECK_EQ(inverse->delta.N0, -expected_delta.N0);
    CHECK_EQ(inverse->delta.N1, -expected_delta.N1);
    CHECK_EQ(inverse->delta.N2, -expected_delta.N2);
    CHECK_EQ(inverse->delta.N3, -expected_delta.N3);
    CHECK_EQ(inverse->delta.N4, -expected_delta.N4);
    return moved->triangulation;
  };

  auto after_two_four = check_expected_delta(
      seed, move_tracker::MoveType4D::TWO_FOUR, S4Counts{0, 1, 4, 5, 2});
  check_expected_delta(after_two_four, move_tracker::MoveType4D::THREE_THREE,
                       S4Counts{});

  check_expected_delta(seed, move_tracker::MoveType4D::TWO_EIGHT,
                       S4Counts{1, 6, 14, 15, 6});

  auto four_six_before = state_with_four_six_site(seed);
  REQUIRE(four_six_before);
  check_expected_delta(*four_six_before, move_tracker::MoveType4D::FOUR_SIX,
                       S4Counts{0, 1, 4, 5, 2});
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

  auto three_three = moves::apply(two_four->triangulation,
                                  move_tracker::MoveType4D::THREE_THREE);
  REQUIRE(three_three);
  check_action_delta(two_four->triangulation, *three_three);

  auto two_eight = moves::apply(seed, move_tracker::MoveType4D::TWO_EIGHT);
  REQUIRE(two_eight);
  check_action_delta(seed, *two_eight);

  auto four_six_before = state_with_four_six_site(seed);
  REQUIRE(four_six_before);
  auto four_six =
      moves::apply(*four_six_before, move_tracker::MoveType4D::FOUR_SIX);
  REQUIRE(four_six);
  check_action_delta(*four_six_before, *four_six);
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
  CHECK_GT(reverse, 0);
}

TEST_CASE("All standard 4D move descriptors expose production proposal sites")
{
  auto triangulation = FoliatedTriangulation4::periodic_seed(4);
  using move_tracker::MoveType4D;
  CHECK_GT(triangulation.candidate_multiplicity(MoveType4D::TWO_FOUR), 0);
  CHECK_GT(triangulation.candidate_multiplicity(MoveType4D::TWO_EIGHT), 0);

  auto two_four = moves::apply(triangulation, MoveType4D::TWO_FOUR);
  REQUIRE(two_four);
  CHECK_GT(two_four->triangulation.candidate_multiplicity(
               MoveType4D::FOUR_TWO),
           0);
  CHECK_GT(two_four->triangulation.candidate_multiplicity(
               MoveType4D::THREE_THREE),
           0);

  auto two_eight = moves::apply(triangulation, MoveType4D::TWO_EIGHT);
  REQUIRE(two_eight);
  CHECK_GT(two_eight->triangulation.candidate_multiplicity(
               MoveType4D::EIGHT_TWO),
           0);

  auto four_six_before = state_with_four_six_site(triangulation);
  REQUIRE(four_six_before);
  CHECK_GT(four_six_before->candidate_multiplicity(MoveType4D::FOUR_SIX), 0);

  auto four_six = moves::apply(*four_six_before, MoveType4D::FOUR_SIX);
  REQUIRE(four_six);
  CHECK_GT(four_six->triangulation.candidate_multiplicity(
               MoveType4D::SIX_FOUR),
           0);
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
