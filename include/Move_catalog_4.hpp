/*******************************************************************************
 Causal Dynamical Triangulations in C++ using CGAL
*******************************************************************************/

/// @file Move_catalog_4.hpp
/// @brief Documented 3+1D CDT move catalogue and proposal accounting.

#ifndef CDT_PLUSPLUS_MOVE_CATALOG_4_HPP
#define CDT_PLUSPLUS_MOVE_CATALOG_4_HPP

#include <array>
#include <optional>
#include <string_view>

#include "Move_tracker.hpp"
#include "S4Action.hpp"

namespace cdt::four_d
{
  enum class ProposalObservable4D
  {
    two_four_sites,
    four_two_sites,
    three_three_sites,
    four_six_sites,
    six_four_sites,
    two_eight_sites,
    eight_two_sites,
    none
  };

  struct ProposalInventory4D
  {
    Int_precision      two_four_sites{0};
    Int_precision      four_two_sites{0};
    Int_precision      three_three_sites{0};
    Int_precision      four_six_sites{0};
    Int_precision      six_four_sites{0};
    Int_precision      two_eight_sites{0};
    Int_precision      eight_two_sites{0};

    [[nodiscard]] auto count(ProposalObservable4D const observable) const
        -> Int_precision
    {
      switch (observable)
      {
        case ProposalObservable4D::two_four_sites: return two_four_sites;
        case ProposalObservable4D::four_two_sites: return four_two_sites;
        case ProposalObservable4D::three_three_sites: return three_three_sites;
        case ProposalObservable4D::four_six_sites: return four_six_sites;
        case ProposalObservable4D::six_four_sites: return six_four_sites;
        case ProposalObservable4D::two_eight_sites: return two_eight_sites;
        case ProposalObservable4D::eight_two_sites: return eight_two_sites;
        case ProposalObservable4D::none: return 0;
      }
      return 0;
    }
  };

  struct MoveDescriptor4D
  {
    move_tracker::MoveType4D move{move_tracker::MoveType4D::NO_MOVE};
    move_tracker::MoveType4D inverse{move_tracker::MoveType4D::NO_MOVE};
    std::string_view         name;
    std::string_view         inverse_name;
    std::string_view         local_subcomplex;
    std::string_view         replacement;
    std::string_view         applicability;
    ProposalObservable4D     proposal_observable{ProposalObservable4D::none};
    S4Counts                 delta;
  };

  [[nodiscard]] constexpr auto all_move_descriptors_4d()
  {
    using enum move_tracker::MoveType4D;
    using enum ProposalObservable4D;
    // clang-format off
    return std::array{
        MoveDescriptor4D{
            .move = TWO_FOUR,
            .inverse = FOUR_TWO,
            .name = "TWO_FOUR",
            .inverse_name = "FOUR_TWO",
            .local_subcomplex =
                "Two adjacent 4-simplices sharing a legal tetrahedral facet.",
            .replacement =
                "Insert the complementary causal edge and replace the pair by "
                "four causal 4-simplices in the same sandwich.",
            .applicability =
                "The shared tetrahedron is internal, causal, and its "
                "replacement keeps all facets paired in periodic time.",
            .proposal_observable = two_four_sites,
            .delta = S4Counts{0, 1, 4, 5, 2, 0, 0, 0, 0, std::nullopt}},
        MoveDescriptor4D{
            .move = FOUR_TWO,
            .inverse = TWO_FOUR,
            .name = "FOUR_TWO",
            .inverse_name = "TWO_FOUR",
            .local_subcomplex =
                "Four 4-simplices around a removable timelike edge.",
            .replacement =
                "Collapse the timelike edge and restore the two-simplex local "
                "subcomplex.",
            .applicability =
                "The edge order is exactly four and the collapse does not "
                "identify distinct boundary vertices.",
            .proposal_observable = four_two_sites,
            .delta = S4Counts{0, -1, -4, -5, -2, 0, 0, 0, 0, std::nullopt}},
        MoveDescriptor4D{
            .move = THREE_THREE,
            .inverse = THREE_THREE,
            .name = "THREE_THREE",
            .inverse_name = "THREE_THREE",
            .local_subcomplex =
                "Three 4-simplices around an order-three timelike triangle.",
            .replacement =
                "Replace the timelike triangle by the complementary timelike "
                "triangle in the same six-vertex ball.",
            .applicability =
                "The dual triangle is absent outside the local ball and the "
                "replacement preserves causal simplex types and gluing.",
            .proposal_observable = three_three_sites,
            .delta = S4Counts{}},
        MoveDescriptor4D{
            .move = FOUR_SIX,
            .inverse = SIX_FOUR,
            .name = "FOUR_SIX",
            .inverse_name = "SIX_FOUR",
            .local_subcomplex =
                "Four 4-simplices around an order-four spatial triangle.",
            .replacement =
                "Replace the spatial triangle by the complementary spatial "
                "edge, producing six 4-simplices.",
            .applicability =
                "The new spatial edge is absent and the four-simplex star has "
                "one adjacent-time apex on each side of the slice.",
            .proposal_observable = four_six_sites,
            .delta = S4Counts{0, 1, 4, 5, 2, 0, 0, 0, 0, std::nullopt}},
        MoveDescriptor4D{
            .move = SIX_FOUR,
            .inverse = FOUR_SIX,
            .name = "SIX_FOUR",
            .inverse_name = "FOUR_SIX",
            .local_subcomplex =
                "Six 4-simplices around an order-six spatial edge.",
            .replacement =
                "Collapse the spatial edge and restore the complementary "
                "order-four spatial triangle.",
            .applicability =
                "The replacement spatial triangle is absent outside the local "
                "ball and the replacement keeps all facets paired.",
            .proposal_observable = six_four_sites,
            .delta = S4Counts{0, -1, -4, -5, -2, 0, 0, 0, 0, std::nullopt}},
        MoveDescriptor4D{
            .move = TWO_EIGHT,
            .inverse = EIGHT_TWO,
            .name = "TWO_EIGHT",
            .inverse_name = "EIGHT_TWO",
            .local_subcomplex =
                "Two 4-simplices sharing a spatial tetrahedron.",
            .replacement =
                "Insert a new spatial vertex in the tetrahedron and replace "
                "the two cones by eight causal 4-simplices.",
            .applicability =
                "The shared tetrahedron has one adjacent-time apex on each "
                "side and the stellar subdivision preserves the closed slice.",
            .proposal_observable = two_eight_sites,
            .delta = S4Counts{1, 6, 14, 15, 6, 0, 0, 0, 0, std::nullopt}},
        MoveDescriptor4D{
            .move = EIGHT_TWO,
            .inverse = TWO_EIGHT,
            .name = "EIGHT_TWO",
            .inverse_name = "TWO_EIGHT",
            .local_subcomplex =
                "Eight 4-simplices around a removable order-eight spatial "
                "vertex.",
            .replacement =
                "Remove the spatial vertex and restore the two causal "
                "4-simplices sharing the spatial tetrahedron.",
            .applicability =
                "The vertex star contains only the expected two cone apices "
                "and the replacement tetrahedron is absent.",
            .proposal_observable = eight_two_sites,
            .delta = S4Counts{-1, -6, -14, -15, -6, 0, 0, 0, 0,
                              std::nullopt}}};
    // clang-format on
  }

  [[nodiscard]] constexpr auto move_descriptor_4d(
      move_tracker::MoveType4D const move) -> MoveDescriptor4D
  {
    for (auto const& descriptor : all_move_descriptors_4d())
    {
      if (descriptor.move == move) { return descriptor; }
    }
    return MoveDescriptor4D{};
  }

  [[nodiscard]] constexpr auto implemented_move_descriptors_4d()
  {
    using enum move_tracker::MoveType4D;
    return std::array{move_descriptor_4d(TWO_FOUR),
                      move_descriptor_4d(FOUR_TWO),
                      move_descriptor_4d(THREE_THREE),
                      move_descriptor_4d(FOUR_SIX),
                      move_descriptor_4d(SIX_FOUR),
                      move_descriptor_4d(TWO_EIGHT),
                      move_descriptor_4d(EIGHT_TWO)};
  }

  [[nodiscard]] constexpr auto reversed_delta(S4Counts const& delta) -> S4Counts
  {
    return S4Counts{-delta.N0,  -delta.N1,   -delta.N2,  -delta.N3,
                    -delta.N4,  -delta.N41,  -delta.N32, -delta.N23,
                    -delta.N14, std::nullopt};
  }

  [[nodiscard]] inline auto proposal_inventory_from_counts(
      S4Counts const& counts) -> ProposalInventory4D
  {
    static_cast<void>(counts);
    return ProposalInventory4D{};
  }

  [[nodiscard]] inline auto local_action_difference(
      S4Counts const& before, S4Counts const& exact_delta,
      S4Couplings const& couplings) -> long double
  {
    auto after = before;
    after.N0 += exact_delta.N0;
    after.N1 += exact_delta.N1;
    after.N2 += exact_delta.N2;
    after.N3 += exact_delta.N3;
    after.N4 += exact_delta.N4;
    after.N41 += exact_delta.N41;
    after.N32 += exact_delta.N32;
    after.N23 += exact_delta.N23;
    after.N14 += exact_delta.N14;
    return S4_action_difference(before, after, couplings);
  }
}  // namespace cdt::four_d

#endif  // CDT_PLUSPLUS_MOVE_CATALOG_4_HPP
