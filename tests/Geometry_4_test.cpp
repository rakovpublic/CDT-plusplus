#include <doctest/doctest.h>

#include "Foliated_triangulation_4.hpp"

using namespace cdt::four_d;
namespace move_tracker = cdt::move_tracker;

TEST_CASE("Persistent 4D seed exposes N0 through N4")
{
  auto triangulation = FoliatedTriangulation4::periodic_seed(3);
  auto counts        = triangulation.counts();

  CHECK_GT(counts.N0, 0);
  CHECK_GT(counts.N1, 0);
  CHECK_GT(counts.N2, 0);
  CHECK_GT(counts.N3, 0);
  CHECK_GT(counts.N4, 0);
}

TEST_CASE("Standard 4D move delta entries preserve Euler characteristic")
{
  auto const delta = FoliatedTriangulation4::move_count_delta(
      move_tracker::MoveType4D::TWO_EIGHT);
  CHECK_EQ(delta.N0, 1);
  CHECK_EQ(delta.N1, 6);
  CHECK_EQ(delta.N2, 14);
  CHECK_EQ(delta.N3, 15);
  CHECK_EQ(delta.N4, 6);
  CHECK_EQ(delta.N0 - delta.N1 + delta.N2 - delta.N3 + delta.N4, 0);
}
