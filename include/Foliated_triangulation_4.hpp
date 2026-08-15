/*******************************************************************************
 Causal Dynamical Triangulations in C++ using CGAL
*******************************************************************************/

/// @file Foliated_triangulation_4.hpp
/// @brief Persistent combinatorial 3+1D CDT triangulation state.

#ifndef CDT_PLUSPLUS_FOLIATED_TRIANGULATION_4_HPP
#define CDT_PLUSPLUS_FOLIATED_TRIANGULATION_4_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Move_catalog_4.hpp"

namespace cdt::four_d
{
  using VertexId  = std::uint64_t;
  using SimplexId = std::uint64_t;

  enum class SimplexType4D
  {
    FOUR_ONE  = 41,
    THREE_TWO = 32,
    TWO_THREE = 23,
    ONE_FOUR  = 14
  };

  struct Vertex4D
  {
    VertexId      id{0};
    Int_precision time{0};
  };

  struct Simplex4D
  {
    SimplexId                               id{0};
    std::array<VertexId, 5>                 vertices{};
    std::array<std::optional<SimplexId>, 5> neighbors{};
    SimplexType4D                           type{SimplexType4D::FOUR_ONE};
  };

  struct ValidationReport
  {
    bool                     standard_cdt_candidate{true};
    std::vector<std::string> errors;

    [[nodiscard]] auto       valid() const -> bool { return errors.empty(); }
  };

  [[nodiscard]] inline auto as_count(SimplexType4D const type) -> Int_precision
  {
    return static_cast<Int_precision>(type);
  }

  /// @brief Persistent combinatorial 3+1D CDT triangulation state.
  ///
  /// The production 4D runner evolves the explicit vertex/simplex incidence
  /// complex. Aggregate counts and spatial profiles are derived caches only;
  /// they are never used as a substitute for local proposal-site enumeration.
  class FoliatedTriangulation4
  {
   public:
    using VertexContainer  = std::vector<Vertex4D>;
    using SimplexContainer = std::vector<Simplex4D>;
    using Profile          = std::vector<Int_precision>;

   private:
    using Edge4D     = std::array<VertexId, 2>;
    using Triangle4D = std::array<VertexId, 3>;
    using Facet4D    = std::array<VertexId, 4>;

    struct TwoFourSite
    {
      Facet4D                  shared_facet{};
      std::array<SimplexId, 2> simplex_ids{};
      Edge4D                   new_edge{};
    };

    struct FourTwoSite
    {
      Edge4D                   edge{};
      std::array<SimplexId, 4> simplex_ids{};
      Facet4D                  replacement_facet{};
    };

    struct ThreeThreeSite
    {
      Triangle4D               triangle{};
      std::array<SimplexId, 3> simplex_ids{};
      Triangle4D               replacement_triangle{};
    };

    struct FourSixSite
    {
      Triangle4D               triangle{};
      std::array<SimplexId, 4> simplex_ids{};
      Edge4D                   new_edge{};
    };

    struct SixFourSite
    {
      Edge4D                   edge{};
      std::array<SimplexId, 6> simplex_ids{};
      Triangle4D               replacement_triangle{};
    };

    struct TwoEightSite
    {
      Facet4D                  shared_facet{};
      std::array<SimplexId, 2> simplex_ids{};
      Int_precision            vertex_time{0};
    };

    struct EightTwoSite
    {
      VertexId                 vertex{0};
      std::array<SimplexId, 8> simplex_ids{};
      Facet4D                  replacement_facet{};
    };

    Int_precision                               m_timeslices{2};
    bool                                        m_periodic{true};
    VertexContainer                             m_vertices;
    SimplexContainer                            m_simplices;
    std::unordered_map<VertexId, Int_precision> m_vertex_times;
    S4Counts                                    m_counts;
    ProposalInventory4D                         m_proposal_inventory;
    Profile                                     m_spatial_profile;
    bool                                        m_closed_s3_slices{true};
    bool                                        m_three_three_forward{true};

    [[nodiscard]] auto vertex_time(VertexId const id) const -> Int_precision
    {
      auto const it = m_vertex_times.find(id);
      return it == m_vertex_times.end() ? -1 : it->second;
    }

    void rebuild_vertex_time_cache()
    {
      m_vertex_times.clear();
      for (auto const& vertex : m_vertices)
      {
        m_vertex_times[vertex.id] = vertex.time;
      }
    }

    [[nodiscard]] auto are_adjacent_times(Int_precision const a,
                                          Int_precision const b) const -> bool
    {
      if (a == b) { return false; }
      if (!m_periodic) { return std::abs(a - b) == 1; }
      auto const delta = (b - a + m_timeslices) % m_timeslices;
      return delta == 1 || delta == m_timeslices - 1;
    }

    [[nodiscard]] auto classify_simplex(std::array<VertexId, 5> const& vertices)
        const -> std::optional<SimplexType4D>
    {
      std::map<Int_precision, int> by_time;
      for (auto const vertex : vertices) { ++by_time[vertex_time(vertex)]; }
      if (by_time.size() != 2) { return std::nullopt; }
      auto first  = by_time.begin();
      auto second = std::next(first);
      if (!are_adjacent_times(first->first, second->first))
      {
        return std::nullopt;
      }

      auto lower_time = first->first;
      auto upper_time = second->first;
      if (m_periodic &&
          ((lower_time + 1) % m_timeslices != upper_time % m_timeslices))
      {
        std::swap(lower_time, upper_time);
      }

      auto lower_vertices = 0;
      for (auto const vertex : vertices)
      {
        if (vertex_time(vertex) == lower_time) { ++lower_vertices; }
      }
      auto const upper_vertices = 5 - lower_vertices;
      if (lower_vertices == 4 && upper_vertices == 1)
      {
        return SimplexType4D::FOUR_ONE;
      }
      if (lower_vertices == 3 && upper_vertices == 2)
      {
        return SimplexType4D::THREE_TWO;
      }
      if (lower_vertices == 2 && upper_vertices == 3)
      {
        return SimplexType4D::TWO_THREE;
      }
      if (lower_vertices == 1 && upper_vertices == 4)
      {
        return SimplexType4D::ONE_FOUR;
      }
      return std::nullopt;
    }

    [[nodiscard]] static auto sorted_vertices(std::array<VertexId, 5> vertices)
    {
      std::ranges::sort(vertices);
      return vertices;
    }

    [[nodiscard]] static auto sorted_edge(Edge4D vertices)
    {
      std::ranges::sort(vertices);
      return vertices;
    }

    [[nodiscard]] static auto sorted_triangle(Triangle4D vertices)
    {
      std::ranges::sort(vertices);
      return vertices;
    }

    [[nodiscard]] static auto sorted_facet(std::array<VertexId, 4> vertices)
    {
      std::ranges::sort(vertices);
      return vertices;
    }

    [[nodiscard]] static auto contains_vertex(
        std::array<VertexId, 5> const& vertices, VertexId const vertex) -> bool
    {
      return std::ranges::find(vertices, vertex) != vertices.end();
    }

    template <std::size_t size>
    [[nodiscard]] static auto contains_vertex(
        std::array<VertexId, size> const& vertices, VertexId const vertex)
        -> bool
    {
      return std::ranges::find(vertices, vertex) != vertices.end();
    }

    [[nodiscard]] auto facet_vertices(Simplex4D const& simplex,
                                      int const omitted_local_index) const
        -> std::array<VertexId, 4>
    {
      std::array<VertexId, 4> facet{};
      auto                    out = 0;
      for (auto index = 0; index < 5; ++index)
      {
        if (index != omitted_local_index)
        {
          facet[static_cast<std::size_t>(out++)] =
              simplex.vertices[static_cast<std::size_t>(index)];
        }
      }
      return sorted_facet(facet);
    }

    [[nodiscard]] auto simplex_by_id(SimplexId const id) const
        -> Simplex4D const*
    {
      auto const it = std::ranges::find_if(
          m_simplices, [id](auto const& simplex) { return simplex.id == id; });
      return it == m_simplices.end() ? nullptr : &*it;
    }

    [[nodiscard]] auto next_simplex_id() const -> SimplexId
    {
      auto next = SimplexId{1};
      for (auto const& simplex : m_simplices)
      {
        next = std::max(next, simplex.id + 1);
      }
      return next;
    }

    [[nodiscard]] auto simplex_contains_edge(Simplex4D const& simplex,
                                             Edge4D const& edge) const -> bool
    {
      return contains_vertex(simplex.vertices, edge[0]) &&
             contains_vertex(simplex.vertices, edge[1]);
    }

    [[nodiscard]] auto simplex_contains_triangle(
        Simplex4D const& simplex, Triangle4D const& triangle) const -> bool
    {
      return contains_vertex(simplex.vertices, triangle[0]) &&
             contains_vertex(simplex.vertices, triangle[1]) &&
             contains_vertex(simplex.vertices, triangle[2]);
    }

    [[nodiscard]] auto edge_exists(Edge4D const& edge) const -> bool
    {
      return std::ranges::any_of(m_simplices, [&](auto const& simplex) {
        return simplex_contains_edge(simplex, edge);
      });
    }

    [[nodiscard]] auto edge_exists_outside(
        Edge4D const& edge, std::set<SimplexId> const& removed_ids) const
        -> bool
    {
      return std::ranges::any_of(m_simplices, [&](auto const& simplex) {
        return !removed_ids.contains(simplex.id) &&
               simplex_contains_edge(simplex, edge);
      });
    }

    [[nodiscard]] auto triangle_exists_outside(
        Triangle4D const& triangle, std::set<SimplexId> const& removed_ids)
        const -> bool
    {
      return std::ranges::any_of(m_simplices, [&](auto const& simplex) {
        return !removed_ids.contains(simplex.id) &&
               simplex_contains_triangle(simplex, triangle);
      });
    }

    template <std::size_t size>
    [[nodiscard]] auto common_time(std::array<VertexId, size> const& vertices)
        const -> std::optional<Int_precision>
    {
      if (vertices.empty()) { return std::nullopt; }
      auto const time = vertex_time(vertices.front());
      if (time < 0) { return std::nullopt; }
      for (auto const vertex : vertices)
      {
        if (vertex_time(vertex) != time) { return std::nullopt; }
      }
      return time;
    }

    [[nodiscard]] auto vertex_order(VertexId const vertex) const
        -> Int_precision
    {
      return static_cast<Int_precision>(std::ranges::count_if(
          m_simplices, [&](auto const& simplex) {
            return contains_vertex(simplex.vertices, vertex);
          }));
    }

    [[nodiscard]] auto simplex_key_exists_outside(
        std::array<VertexId, 5>    vertices,
        std::set<SimplexId> const& removed_ids) const -> bool
    {
      auto const key = sorted_vertices(vertices);
      return std::ranges::any_of(m_simplices, [&](auto const& simplex) {
        return !removed_ids.contains(simplex.id) &&
               sorted_vertices(simplex.vertices) == key;
      });
    }

    [[nodiscard]] auto recompute_counts_from_complex() const -> S4Counts
    {
      S4Counts counts;
      counts.N0 = static_cast<Int_precision>(m_vertices.size());
      counts.N4 = static_cast<Int_precision>(m_simplices.size());

      std::set<std::array<VertexId, 2>> edges;
      std::set<std::array<VertexId, 2>> timelike_edges;
      std::set<std::array<VertexId, 3>> triangles;
      std::set<std::array<VertexId, 3>> mixed_triangles;
      std::set<std::array<VertexId, 4>> tetrahedra;
      std::set<std::array<VertexId, 4>> spatial_tetrahedra;
      std::set<std::array<VertexId, 4>> timelike_tetrahedra;

      for (auto const& simplex : m_simplices)
      {
        switch (simplex.type)
        {
          case SimplexType4D::FOUR_ONE: ++counts.N41; break;
          case SimplexType4D::THREE_TWO: ++counts.N32; break;
          case SimplexType4D::TWO_THREE: ++counts.N23; break;
          case SimplexType4D::ONE_FOUR: ++counts.N14; break;
        }

        for (auto i = 0; i < 5; ++i)
        {
          for (auto j = i + 1; j < 5; ++j)
          {
            auto edge = std::array{simplex.vertices[static_cast<size_t>(i)],
                                   simplex.vertices[static_cast<size_t>(j)]};
            std::ranges::sort(edge);
            edges.insert(edge);
            if (vertex_time(edge[0]) != vertex_time(edge[1]))
            {
              timelike_edges.insert(edge);
            }
          }
        }
        for (auto i = 0; i < 5; ++i)
        {
          for (auto j = i + 1; j < 5; ++j)
          {
            for (auto k = j + 1; k < 5; ++k)
            {
              auto triangle =
                  std::array{simplex.vertices[static_cast<size_t>(i)],
                             simplex.vertices[static_cast<size_t>(j)],
                             simplex.vertices[static_cast<size_t>(k)]};
              std::ranges::sort(triangle);
              triangles.insert(triangle);
              std::set<Int_precision> times;
              for (auto const vertex : triangle)
              {
                times.insert(vertex_time(vertex));
              }
              if (times.size() > 1) { mixed_triangles.insert(triangle); }
            }
          }
        }
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          auto const facet = facet_vertices(simplex, omitted);
          tetrahedra.insert(facet);
          std::set<Int_precision> times;
          for (auto const vertex : facet) { times.insert(vertex_time(vertex)); }
          if (times.size() == 1) { spatial_tetrahedra.insert(facet); }
          else { timelike_tetrahedra.insert(facet); }
        }
      }

      counts.N1             = static_cast<Int_precision>(edges.size());
      counts.N2             = static_cast<Int_precision>(triangles.size());
      counts.N3             = static_cast<Int_precision>(tetrahedra.size());
      counts.class_resolved = S4ClassResolvedCounts{
          static_cast<Int_precision>(spatial_tetrahedra.size()),
          static_cast<Int_precision>(timelike_edges.size()),
          static_cast<Int_precision>(mixed_triangles.size()),
          static_cast<Int_precision>(timelike_tetrahedra.size())};
      return counts;
    }

    [[nodiscard]] auto recompute_spatial_profile() const -> Profile
    {
      Profile profile(static_cast<std::size_t>(m_timeslices), 0);
      std::set<std::array<VertexId, 4>> seen;
      for (auto const& simplex : m_simplices)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          auto const facet = facet_vertices(simplex, omitted);
          if (seen.contains(facet)) { continue; }
          seen.insert(facet);
          std::set<Int_precision> times;
          for (auto const vertex : facet) { times.insert(vertex_time(vertex)); }
          if (times.size() == 1)
          {
            auto const time = *times.begin();
            if (time >= 0 && time < m_timeslices)
            {
              ++profile[static_cast<std::size_t>(time)];
            }
          }
        }
      }
      return profile;
    }

    void rebuild_neighbors()
    {
      for (auto& simplex : m_simplices)
      {
        simplex.neighbors.fill(std::nullopt);
      }
      std::map<Facet4D, std::vector<std::pair<std::size_t, int>>> incidence;
      for (std::size_t simplex_index = 0; simplex_index < m_simplices.size();
           ++simplex_index)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          incidence[facet_vertices(m_simplices[simplex_index], omitted)]
              .push_back({simplex_index, omitted});
        }
      }
      for (auto const& [_, incident] : incidence)
      {
        if (incident.size() != 2) { continue; }
        auto const [first_index, first_local]   = incident[0];
        auto const [second_index, second_local] = incident[1];
        m_simplices[first_index]
            .neighbors[static_cast<std::size_t>(first_local)] =
            m_simplices[second_index].id;
        m_simplices[second_index]
            .neighbors[static_cast<std::size_t>(second_local)] =
            m_simplices[first_index].id;
      }
    }

    [[nodiscard]] auto make_site_inventory() const -> ProposalInventory4D
    {
      return ProposalInventory4D{
          static_cast<Int_precision>(enumerate_two_four_sites().size()),
          static_cast<Int_precision>(enumerate_four_two_sites().size()),
          static_cast<Int_precision>(enumerate_three_three_sites().size()),
          static_cast<Int_precision>(enumerate_four_six_sites().size()),
          static_cast<Int_precision>(enumerate_six_four_sites().size()),
          static_cast<Int_precision>(enumerate_two_eight_sites().size()),
          static_cast<Int_precision>(enumerate_eight_two_sites().size())};
    }

    void refresh_derived_state()
    {
      rebuild_vertex_time_cache();
      rebuild_neighbors();
      m_counts             = recompute_counts_from_complex();
      m_spatial_profile    = recompute_spatial_profile();
      m_proposal_inventory = make_site_inventory();
    }

    [[nodiscard]] auto replacement_is_legal(
        std::vector<std::array<VertexId, 5>> const& replacement,
        std::set<SimplexId> const&                  removed_ids) const -> bool
    {
      std::set<std::array<VertexId, 5>> replacement_keys;
      for (auto const& vertices : replacement)
      {
        if (std::set<VertexId>(vertices.begin(), vertices.end()).size() != 5)
        {
          return false;
        }
        if (!classify_simplex(vertices)) { return false; }
        auto const key = sorted_vertices(vertices);
        if (!replacement_keys.insert(key).second) { return false; }
        if (simplex_key_exists_outside(vertices, removed_ids)) { return false; }
      }
      return true;
    }

    [[nodiscard]] static auto cluster_boundary_facets(
        std::vector<std::array<VertexId, 5>> const& cluster)
        -> std::optional<std::set<Facet4D>>
    {
      std::map<Facet4D, int> facet_incidence;
      for (auto const& simplex : cluster)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          Facet4D facet{};
          auto    out = 0;
          for (auto index = 0; index < 5; ++index)
          {
            if (index == omitted) { continue; }
            facet[static_cast<std::size_t>(out++)] =
                simplex[static_cast<std::size_t>(index)];
          }
          ++facet_incidence[sorted_facet(facet)];
        }
      }

      std::set<Facet4D> boundary;
      for (auto const& [facet, incidence] : facet_incidence)
      {
        if (incidence == 1) { boundary.insert(facet); }
        else if (incidence != 2) { return std::nullopt; }
      }
      return boundary;
    }

    [[nodiscard]] auto replacement_preserves_valid_complex(
        std::set<SimplexId> const&                  removed_ids,
        std::vector<std::array<VertexId, 5>> const& replacement) const -> bool
    {
      if (!replacement_is_legal(replacement, removed_ids)) { return false; }

      std::vector<std::array<VertexId, 5>> removed_cluster;
      removed_cluster.reserve(removed_ids.size());
      for (auto const& simplex : m_simplices)
      {
        if (removed_ids.contains(simplex.id))
        {
          removed_cluster.push_back(simplex.vertices);
        }
      }
      if (removed_cluster.size() != removed_ids.size()) { return false; }

      auto const removed_boundary = cluster_boundary_facets(removed_cluster);
      auto const replacement_boundary = cluster_boundary_facets(replacement);
      if (!removed_boundary || !replacement_boundary ||
          *removed_boundary != *replacement_boundary)
      {
        return false;
      }

      std::map<Facet4D, int> outside_incidence;
      for (auto const& simplex : m_simplices)
      {
        if (removed_ids.contains(simplex.id)) { continue; }
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          ++outside_incidence[facet_vertices(simplex, omitted)];
        }
      }

      std::map<Facet4D, int> replacement_incidence;
      for (auto const& simplex : replacement)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          Facet4D facet{};
          auto    out = 0;
          for (auto index = 0; index < 5; ++index)
          {
            if (index == omitted) { continue; }
            facet[static_cast<std::size_t>(out++)] =
                simplex[static_cast<std::size_t>(index)];
          }
          ++replacement_incidence[sorted_facet(facet)];
        }
      }

      for (auto const& [facet, incidence] : replacement_incidence)
      {
        auto const outside_count =
            outside_incidence.contains(facet) ? outside_incidence[facet] : 0;
        if (incidence == 2 && outside_count != 0) { return false; }
        if (incidence == 1 && outside_count != 1) { return false; }
      }
      return true;
    }

    [[nodiscard]] auto enumerate_two_four_sites() const
        -> std::vector<TwoFourSite>
    {
      std::map<Facet4D, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          incidence[facet_vertices(simplex, omitted)].push_back(simplex.id);
        }
      }

      std::vector<TwoFourSite> sites;
      for (auto const& [facet, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 2) { continue; }
        auto const* first  = simplex_by_id(simplex_ids[0]);
        auto const* second = simplex_by_id(simplex_ids[1]);
        if (first == nullptr || second == nullptr) { continue; }

        std::optional<VertexId> first_apex;
        std::optional<VertexId> second_apex;
        for (auto const vertex : first->vertices)
        {
          if (!contains_vertex(facet, vertex)) { first_apex = vertex; }
        }
        for (auto const vertex : second->vertices)
        {
          if (!contains_vertex(facet, vertex)) { second_apex = vertex; }
        }
        if (!first_apex || !second_apex || *first_apex == *second_apex)
        {
          continue;
        }

        auto const new_edge = sorted_edge(Edge4D{*first_apex, *second_apex});
        if (edge_exists(new_edge)) { continue; }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(4);
        for (auto omitted_vertex : facet)
        {
          std::array<VertexId, 5> vertices{};
          auto                    out = 0;
          for (auto const vertex : facet)
          {
            if (vertex != omitted_vertex)
            {
              vertices[static_cast<std::size_t>(out++)] = vertex;
            }
          }
          vertices[static_cast<std::size_t>(out++)] = *first_apex;
          vertices[static_cast<std::size_t>(out++)] = *second_apex;
          replacement.push_back(vertices);
        }
        auto const removed_ids =
            std::set<SimplexId>{simplex_ids[0], simplex_ids[1]};
        if (!replacement_preserves_valid_complex(removed_ids, replacement))
        {
          continue;
        }

        auto site         = TwoFourSite{};
        site.shared_facet = facet;
        site.simplex_ids  = std::array{simplex_ids[0], simplex_ids[1]};
        site.new_edge     = new_edge;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.shared_facet, site.new_edge, site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto enumerate_four_two_sites() const
        -> std::vector<FourTwoSite>
    {
      std::map<Edge4D, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto i = 0; i < 5; ++i)
        {
          for (auto j = i + 1; j < 5; ++j)
          {
            incidence[sorted_edge(Edge4D{
                          simplex.vertices[static_cast<std::size_t>(i)],
                          simplex.vertices[static_cast<std::size_t>(j)]})]
                .push_back(simplex.id);
          }
        }
      }

      std::vector<FourTwoSite> sites;
      for (auto const& [edge, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 4) { continue; }
        if (vertex_time(edge[0]) == vertex_time(edge[1]) ||
            !are_adjacent_times(vertex_time(edge[0]), vertex_time(edge[1])))
        {
          continue;
        }
        std::set<VertexId> union_vertices;
        auto               all_simplices_found = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr || !simplex_contains_edge(*simplex, edge))
          {
            all_simplices_found = false;
            break;
          }
          union_vertices.insert(simplex->vertices.begin(),
                                simplex->vertices.end());
        }
        if (!all_simplices_found || union_vertices.size() != 6) { continue; }

        Facet4D replacement_facet{};
        auto    out = 0;
        for (auto const vertex : union_vertices)
        {
          if (vertex != edge[0] && vertex != edge[1])
          {
            replacement_facet[static_cast<std::size_t>(out++)] = vertex;
          }
        }
        if (out != 4) { continue; }
        replacement_facet = sorted_facet(replacement_facet);

        std::set<VertexId> omitted_vertices;
        auto               star_has_expected_form = true;
        for (auto const id : simplex_ids)
        {
          auto const*             simplex                  = simplex_by_id(id);
          auto                    link_vertices_in_simplex = 0;
          std::optional<VertexId> omitted;
          for (auto const vertex : replacement_facet)
          {
            if (contains_vertex(simplex->vertices, vertex))
            {
              ++link_vertices_in_simplex;
            }
            else { omitted = vertex; }
          }
          if (link_vertices_in_simplex != 3 || !omitted)
          {
            star_has_expected_form = false;
            break;
          }
          omitted_vertices.insert(*omitted);
        }
        if (!star_has_expected_form || omitted_vertices.size() != 4)
        {
          continue;
        }

        auto first_replacement = std::array<VertexId, 5>{
            replacement_facet[0], replacement_facet[1], replacement_facet[2],
            replacement_facet[3], edge[0]};
        auto second_replacement = std::array<VertexId, 5>{
            replacement_facet[0], replacement_facet[1], replacement_facet[2],
            replacement_facet[3], edge[1]};
        auto const removed_ids =
            std::set<SimplexId>(simplex_ids.begin(), simplex_ids.end());
        auto const replacements =
            std::vector{first_replacement, second_replacement};
        if (!replacement_preserves_valid_complex(removed_ids, replacements))
        {
          continue;
        }

        auto site              = FourTwoSite{};
        site.edge              = edge;
        site.simplex_ids       = std::array{simplex_ids[0], simplex_ids[1],
                                      simplex_ids[2], simplex_ids[3]};
        site.replacement_facet = replacement_facet;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.edge, site.replacement_facet, site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto enumerate_three_three_sites() const
        -> std::vector<ThreeThreeSite>
    {
      std::map<Triangle4D, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto i = 0; i < 5; ++i)
        {
          for (auto j = i + 1; j < 5; ++j)
          {
            for (auto k = j + 1; k < 5; ++k)
            {
              incidence[sorted_triangle(Triangle4D{
                            simplex.vertices[static_cast<std::size_t>(i)],
                            simplex.vertices[static_cast<std::size_t>(j)],
                            simplex.vertices[static_cast<std::size_t>(k)]})]
                  .push_back(simplex.id);
            }
          }
        }
      }

      std::vector<ThreeThreeSite> sites;
      for (auto const& [triangle, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 3) { continue; }
        if (common_time(triangle)) { continue; }

        std::set<VertexId> union_vertices;
        auto               star_has_triangle = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr ||
              !simplex_contains_triangle(*simplex, triangle))
          {
            star_has_triangle = false;
            break;
          }
          union_vertices.insert(simplex->vertices.begin(),
                                simplex->vertices.end());
        }
        if (!star_has_triangle || union_vertices.size() != 6) { continue; }

        Triangle4D replacement_triangle{};
        auto       out = 0;
        for (auto const vertex : union_vertices)
        {
          if (!contains_vertex(triangle, vertex))
          {
            replacement_triangle[static_cast<std::size_t>(out++)] = vertex;
          }
        }
        if (out != 3 || common_time(replacement_triangle)) { continue; }
        replacement_triangle = sorted_triangle(replacement_triangle);

        auto const removed_ids =
            std::set<SimplexId>(simplex_ids.begin(), simplex_ids.end());
        if (triangle_exists_outside(replacement_triangle, removed_ids))
        {
          continue;
        }

        std::set<VertexId> omitted_dual_vertices;
        auto               star_has_expected_form = true;
        for (auto const id : simplex_ids)
        {
          auto const*             simplex = simplex_by_id(id);
          auto                    dual_vertices_in_simplex = 0;
          std::optional<VertexId> omitted;
          for (auto const vertex : replacement_triangle)
          {
            if (contains_vertex(simplex->vertices, vertex))
            {
              ++dual_vertices_in_simplex;
            }
            else { omitted = vertex; }
          }
          if (dual_vertices_in_simplex != 2 || !omitted)
          {
            star_has_expected_form = false;
            break;
          }
          omitted_dual_vertices.insert(*omitted);
        }
        if (!star_has_expected_form || omitted_dual_vertices.size() != 3)
        {
          continue;
        }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(3);
        for (auto omitted_vertex : triangle)
        {
          std::array<VertexId, 5> vertices{};
          auto                    vertex_out = 0;
          for (auto const vertex : replacement_triangle)
          {
            vertices[static_cast<std::size_t>(vertex_out++)] = vertex;
          }
          for (auto const vertex : triangle)
          {
            if (vertex != omitted_vertex)
            {
              vertices[static_cast<std::size_t>(vertex_out++)] = vertex;
            }
          }
          replacement.push_back(vertices);
        }
        if (!replacement_preserves_valid_complex(removed_ids, replacement))
        {
          continue;
        }

        auto site                  = ThreeThreeSite{};
        site.triangle              = triangle;
        site.simplex_ids           = std::array{simplex_ids[0], simplex_ids[1],
                                      simplex_ids[2]};
        site.replacement_triangle  = replacement_triangle;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.triangle, site.replacement_triangle,
                          site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto enumerate_four_six_sites() const
        -> std::vector<FourSixSite>
    {
      std::map<Triangle4D, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto i = 0; i < 5; ++i)
        {
          for (auto j = i + 1; j < 5; ++j)
          {
            for (auto k = j + 1; k < 5; ++k)
            {
              incidence[sorted_triangle(Triangle4D{
                            simplex.vertices[static_cast<std::size_t>(i)],
                            simplex.vertices[static_cast<std::size_t>(j)],
                            simplex.vertices[static_cast<std::size_t>(k)]})]
                  .push_back(simplex.id);
            }
          }
        }
      }

      std::vector<FourSixSite> sites;
      for (auto const& [triangle, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 4) { continue; }
        auto const slice_time = common_time(triangle);
        if (!slice_time) { continue; }

        std::set<VertexId> union_vertices;
        auto               star_has_triangle = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr ||
              !simplex_contains_triangle(*simplex, triangle))
          {
            star_has_triangle = false;
            break;
          }
          union_vertices.insert(simplex->vertices.begin(),
                                simplex->vertices.end());
        }
        if (!star_has_triangle || union_vertices.size() != 7) { continue; }

        std::vector<VertexId> spatial_edge_vertices;
        std::vector<VertexId> apices;
        for (auto const vertex : union_vertices)
        {
          if (contains_vertex(triangle, vertex)) { continue; }
          auto const time = vertex_time(vertex);
          if (time == *slice_time) { spatial_edge_vertices.push_back(vertex); }
          else if (are_adjacent_times(time, *slice_time))
          {
            apices.push_back(vertex);
          }
        }
        if (spatial_edge_vertices.size() != 2 || apices.size() != 2)
        {
          continue;
        }
        if (m_timeslices > 2 &&
            vertex_time(apices[0]) == vertex_time(apices[1]))
        {
          continue;
        }

        auto const new_edge =
            sorted_edge(Edge4D{spatial_edge_vertices[0],
                               spatial_edge_vertices[1]});
        auto const removed_ids =
            std::set<SimplexId>(simplex_ids.begin(), simplex_ids.end());
        if (edge_exists_outside(new_edge, removed_ids)) { continue; }

        std::set<std::pair<VertexId, VertexId>> edge_apex_pairs;
        auto                                    star_has_expected_form = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          std::vector<VertexId> edge_vertices_in_simplex;
          std::vector<VertexId> apices_in_simplex;
          for (auto const vertex : spatial_edge_vertices)
          {
            if (contains_vertex(simplex->vertices, vertex))
            {
              edge_vertices_in_simplex.push_back(vertex);
            }
          }
          for (auto const vertex : apices)
          {
            if (contains_vertex(simplex->vertices, vertex))
            {
              apices_in_simplex.push_back(vertex);
            }
          }
          if (edge_vertices_in_simplex.size() != 1 ||
              apices_in_simplex.size() != 1)
          {
            star_has_expected_form = false;
            break;
          }
          edge_apex_pairs.insert(
              {edge_vertices_in_simplex.front(), apices_in_simplex.front()});
        }
        if (!star_has_expected_form || edge_apex_pairs.size() != 4)
        {
          continue;
        }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(6);
        for (auto const apex : apices)
        {
          for (auto const omitted_vertex : triangle)
          {
            std::array<VertexId, 5> vertices{};
            auto                    out = 0;
            vertices[static_cast<std::size_t>(out++)] = new_edge[0];
            vertices[static_cast<std::size_t>(out++)] = new_edge[1];
            vertices[static_cast<std::size_t>(out++)] = apex;
            for (auto const vertex : triangle)
            {
              if (vertex != omitted_vertex)
              {
                vertices[static_cast<std::size_t>(out++)] = vertex;
              }
            }
            replacement.push_back(vertices);
          }
        }
        if (!replacement_preserves_valid_complex(removed_ids, replacement))
        {
          continue;
        }

        auto site         = FourSixSite{};
        site.triangle     = triangle;
        site.simplex_ids  = std::array{simplex_ids[0], simplex_ids[1],
                                      simplex_ids[2], simplex_ids[3]};
        site.new_edge     = new_edge;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.triangle, site.new_edge, site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto enumerate_six_four_sites() const
        -> std::vector<SixFourSite>
    {
      std::map<Edge4D, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto i = 0; i < 5; ++i)
        {
          for (auto j = i + 1; j < 5; ++j)
          {
            incidence[sorted_edge(Edge4D{
                          simplex.vertices[static_cast<std::size_t>(i)],
                          simplex.vertices[static_cast<std::size_t>(j)]})]
                .push_back(simplex.id);
          }
        }
      }

      std::vector<SixFourSite> sites;
      for (auto const& [edge, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 6) { continue; }
        auto const slice_time = common_time(edge);
        if (!slice_time) { continue; }

        std::set<VertexId> union_vertices;
        auto               star_has_edge = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr || !simplex_contains_edge(*simplex, edge))
          {
            star_has_edge = false;
            break;
          }
          union_vertices.insert(simplex->vertices.begin(),
                                simplex->vertices.end());
        }
        if (!star_has_edge || union_vertices.size() != 7) { continue; }

        std::vector<VertexId> triangle_vertices;
        std::vector<VertexId> apices;
        for (auto const vertex : union_vertices)
        {
          if (contains_vertex(edge, vertex)) { continue; }
          auto const time = vertex_time(vertex);
          if (time == *slice_time) { triangle_vertices.push_back(vertex); }
          else if (are_adjacent_times(time, *slice_time))
          {
            apices.push_back(vertex);
          }
        }
        if (triangle_vertices.size() != 3 || apices.size() != 2)
        {
          continue;
        }
        if (m_timeslices > 2 &&
            vertex_time(apices[0]) == vertex_time(apices[1]))
        {
          continue;
        }

        auto const replacement_triangle =
            sorted_triangle(Triangle4D{triangle_vertices[0],
                                       triangle_vertices[1],
                                       triangle_vertices[2]});
        auto const removed_ids =
            std::set<SimplexId>(simplex_ids.begin(), simplex_ids.end());
        if (triangle_exists_outside(replacement_triangle, removed_ids))
        {
          continue;
        }

        std::set<std::pair<VertexId, VertexId>> omitted_apex_pairs;
        auto                                    star_has_expected_form = true;
        for (auto const id : simplex_ids)
        {
          auto const*             simplex = simplex_by_id(id);
          std::optional<VertexId> omitted_triangle_vertex;
          std::optional<VertexId> apex_in_simplex;
          auto                    triangle_count = 0;
          for (auto const vertex : replacement_triangle)
          {
            if (contains_vertex(simplex->vertices, vertex)) { ++triangle_count; }
            else { omitted_triangle_vertex = vertex; }
          }
          for (auto const vertex : apices)
          {
            if (contains_vertex(simplex->vertices, vertex))
            {
              apex_in_simplex = vertex;
            }
          }
          if (triangle_count != 2 || !omitted_triangle_vertex ||
              !apex_in_simplex)
          {
            star_has_expected_form = false;
            break;
          }
          omitted_apex_pairs.insert(
              {*omitted_triangle_vertex, *apex_in_simplex});
        }
        if (!star_has_expected_form || omitted_apex_pairs.size() != 6)
        {
          continue;
        }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(4);
        for (auto const endpoint : edge)
        {
          for (auto const apex : apices)
          {
            replacement.push_back(std::array<VertexId, 5>{
                replacement_triangle[0], replacement_triangle[1],
                replacement_triangle[2], endpoint, apex});
          }
        }
        if (!replacement_preserves_valid_complex(removed_ids, replacement))
        {
          continue;
        }

        auto site                  = SixFourSite{};
        site.edge                  = edge;
        site.simplex_ids           = std::array{simplex_ids[0], simplex_ids[1],
                                      simplex_ids[2], simplex_ids[3],
                                      simplex_ids[4], simplex_ids[5]};
        site.replacement_triangle  = replacement_triangle;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.edge, site.replacement_triangle,
                          site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto enumerate_two_eight_sites() const
        -> std::vector<TwoEightSite>
    {
      std::map<Facet4D, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          incidence[facet_vertices(simplex, omitted)].push_back(simplex.id);
        }
      }

      std::vector<TwoEightSite> sites;
      for (auto const& [facet, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 2) { continue; }
        auto const slice_time = common_time(facet);
        if (!slice_time) { continue; }

        std::vector<VertexId> apices;
        auto                  star_has_facet = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr)
          {
            star_has_facet = false;
            break;
          }
          std::optional<VertexId> apex;
          for (auto const vertex : simplex->vertices)
          {
            if (!contains_vertex(facet, vertex)) { apex = vertex; }
          }
          if (!apex || !are_adjacent_times(vertex_time(*apex), *slice_time))
          {
            star_has_facet = false;
            break;
          }
          apices.push_back(*apex);
        }
        if (!star_has_facet || apices.size() != 2 || apices[0] == apices[1])
        {
          continue;
        }
        if (m_timeslices > 2 &&
            vertex_time(apices[0]) == vertex_time(apices[1]))
        {
          continue;
        }

        auto const new_vertex = next_vertex_id();
        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(8);
        for (auto const apex : apices)
        {
          for (auto const omitted_vertex : facet)
          {
            std::array<VertexId, 5> vertices{};
            auto                    out = 0;
            vertices[static_cast<std::size_t>(out++)] = new_vertex;
            vertices[static_cast<std::size_t>(out++)] = apex;
            for (auto const vertex : facet)
            {
              if (vertex != omitted_vertex)
              {
                vertices[static_cast<std::size_t>(out++)] = vertex;
              }
            }
            replacement.push_back(vertices);
          }
        }

        auto probe = *this;
        probe.m_vertices.push_back(Vertex4D{new_vertex, *slice_time});
        probe.rebuild_vertex_time_cache();
        auto const removed_ids =
            std::set<SimplexId>{simplex_ids[0], simplex_ids[1]};
        if (!probe.replacement_preserves_valid_complex(removed_ids,
                                                       replacement))
        {
          continue;
        }

        auto site         = TwoEightSite{};
        site.shared_facet = facet;
        site.simplex_ids  = std::array{simplex_ids[0], simplex_ids[1]};
        site.vertex_time  = *slice_time;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.shared_facet, site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto enumerate_eight_two_sites() const
        -> std::vector<EightTwoSite>
    {
      std::map<VertexId, std::vector<SimplexId>> incidence;
      for (auto const& simplex : m_simplices)
      {
        for (auto const vertex : simplex.vertices)
        {
          incidence[vertex].push_back(simplex.id);
        }
      }

      std::vector<EightTwoSite> sites;
      for (auto const& [vertex, simplex_ids] : incidence)
      {
        if (simplex_ids.size() != 8 || vertex_order(vertex) != 8)
        {
          continue;
        }
        auto const slice_time = vertex_time(vertex);
        if (slice_time < 0) { continue; }

        std::set<VertexId> union_vertices;
        auto               star_has_vertex = true;
        for (auto const id : simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr || !contains_vertex(simplex->vertices, vertex))
          {
            star_has_vertex = false;
            break;
          }
          union_vertices.insert(simplex->vertices.begin(),
                                simplex->vertices.end());
        }
        if (!star_has_vertex || union_vertices.size() != 7) { continue; }

        std::vector<VertexId> facet_vertices;
        std::vector<VertexId> apices;
        for (auto const candidate : union_vertices)
        {
          if (candidate == vertex) { continue; }
          auto const time = vertex_time(candidate);
          if (time == slice_time) { facet_vertices.push_back(candidate); }
          else if (are_adjacent_times(time, slice_time))
          {
            apices.push_back(candidate);
          }
        }
        if (facet_vertices.size() != 4 || apices.size() != 2)
        {
          continue;
        }
        if (m_timeslices > 2 &&
            vertex_time(apices[0]) == vertex_time(apices[1]))
        {
          continue;
        }

        auto replacement_facet =
            sorted_facet(Facet4D{facet_vertices[0], facet_vertices[1],
                                 facet_vertices[2], facet_vertices[3]});
        std::set<std::pair<VertexId, VertexId>> omitted_apex_pairs;
        auto                                    star_has_expected_form = true;
        for (auto const id : simplex_ids)
        {
          auto const*             simplex = simplex_by_id(id);
          std::optional<VertexId> omitted_facet_vertex;
          std::optional<VertexId> apex_in_simplex;
          auto                    facet_count = 0;
          for (auto const facet_vertex : replacement_facet)
          {
            if (contains_vertex(simplex->vertices, facet_vertex))
            {
              ++facet_count;
            }
            else { omitted_facet_vertex = facet_vertex; }
          }
          for (auto const apex : apices)
          {
            if (contains_vertex(simplex->vertices, apex))
            {
              apex_in_simplex = apex;
            }
          }
          if (facet_count != 3 || !omitted_facet_vertex || !apex_in_simplex)
          {
            star_has_expected_form = false;
            break;
          }
          omitted_apex_pairs.insert({*omitted_facet_vertex, *apex_in_simplex});
        }
        if (!star_has_expected_form || omitted_apex_pairs.size() != 8)
        {
          continue;
        }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(2);
        for (auto const apex : apices)
        {
          replacement.push_back(std::array<VertexId, 5>{
              replacement_facet[0], replacement_facet[1],
              replacement_facet[2], replacement_facet[3], apex});
        }
        auto const removed_ids =
            std::set<SimplexId>(simplex_ids.begin(), simplex_ids.end());
        if (!replacement_preserves_valid_complex(removed_ids, replacement))
        {
          continue;
        }

        auto site              = EightTwoSite{};
        site.vertex            = vertex;
        site.simplex_ids       = std::array{
            simplex_ids[0], simplex_ids[1], simplex_ids[2], simplex_ids[3],
            simplex_ids[4], simplex_ids[5], simplex_ids[6], simplex_ids[7]};
        site.replacement_facet = replacement_facet;
        sites.push_back(site);
      }
      std::ranges::sort(sites, {}, [](auto const& site) {
        return std::tuple{site.vertex, site.replacement_facet,
                          site.simplex_ids};
      });
      return sites;
    }

    [[nodiscard]] auto replace_simplices_and_vertices(
        std::set<SimplexId> const&                  removed_ids,
        std::vector<std::array<VertexId, 5>> const& replacement,
        std::vector<Vertex4D> const&                added_vertices = {},
        std::set<VertexId> const&                   removed_vertices = {})
        -> bool
    {
      auto before = *this;
      for (auto const& vertex : added_vertices)
      {
        if (vertex.time < 0 || vertex.time >= m_timeslices ||
            m_vertex_times.contains(vertex.id))
        {
          *this = before;
          return false;
        }
        m_vertices.push_back(vertex);
        m_vertex_times[vertex.id] = vertex.time;
      }

      if (!replacement_preserves_valid_complex(removed_ids, replacement))
      {
        *this = before;
        return false;
      }

      SimplexContainer updated;
      updated.reserve(m_simplices.size() - removed_ids.size() +
                      replacement.size());
      for (auto const& simplex : m_simplices)
      {
        if (removed_ids.contains(simplex.id)) { continue; }
        if (std::ranges::any_of(simplex.vertices, [&](auto const vertex) {
              return removed_vertices.contains(vertex);
            }))
        {
          *this = before;
          return false;
        }
        updated.push_back(simplex);
      }

      if (!removed_vertices.empty())
      {
        auto const erase_begin = std::ranges::remove_if(
            m_vertices, [&](auto const& vertex) {
              return removed_vertices.contains(vertex.id);
            });
        m_vertices.erase(erase_begin.begin(), erase_begin.end());
        for (auto const vertex : removed_vertices)
        {
          m_vertex_times.erase(vertex);
        }
      }

      auto next_id = next_simplex_id();
      for (auto const& vertices : replacement)
      {
        auto type = classify_simplex(vertices);
        if (!type)
        {
          *this = before;
          return false;
        }
        Simplex4D simplex;
        simplex.id       = next_id++;
        simplex.vertices = vertices;
        simplex.type     = *type;
        updated.push_back(simplex);
      }

      m_simplices = std::move(updated);
      refresh_derived_state();
      if (!validate().valid())
      {
        *this = before;
        return false;
      }
      return true;
    }

    [[nodiscard]] auto replace_simplices(
        std::set<SimplexId> const&                  removed_ids,
        std::vector<std::array<VertexId, 5>> const& replacement) -> bool
    {
      return replace_simplices_and_vertices(removed_ids, replacement);
    }

    [[nodiscard]] auto spacelike_facets_by_slice() const
        -> std::vector<std::set<std::array<VertexId, 4>>>
    {
      std::vector<std::set<std::array<VertexId, 4>>> facets(
          static_cast<std::size_t>(m_timeslices));
      std::set<std::array<VertexId, 4>> seen;
      for (auto const& simplex : m_simplices)
      {
        for (auto omitted = 0; omitted < 5; ++omitted)
        {
          auto const facet = facet_vertices(simplex, omitted);
          if (!seen.insert(facet).second) { continue; }
          std::set<Int_precision> times;
          for (auto const vertex : facet) { times.insert(vertex_time(vertex)); }
          if (times.size() != 1) { continue; }
          auto const time = *times.begin();
          if (time >= 0 && time < m_timeslices)
          {
            facets[static_cast<std::size_t>(time)].insert(facet);
          }
        }
      }
      return facets;
    }

    [[nodiscard]] auto derived_slice_euler_characteristics() const -> Profile
    {
      if (m_simplices.empty())
      {
        return Profile(static_cast<std::size_t>(m_timeslices),
                       static_cast<Int_precision>(1));
      }

      Profile    result(static_cast<std::size_t>(m_timeslices), 0);
      auto const facets_by_slice = spacelike_facets_by_slice();
      for (std::size_t slice = 0; slice < facets_by_slice.size(); ++slice)
      {
        std::set<VertexId>                vertices;
        std::set<std::array<VertexId, 2>> edges;
        std::set<std::array<VertexId, 3>> triangles;
        for (auto const& facet : facets_by_slice[slice])
        {
          for (auto const vertex : facet) { vertices.insert(vertex); }
          for (auto i = 0; i < 4; ++i)
          {
            for (auto j = i + 1; j < 4; ++j)
            {
              auto edge = std::array{facet[static_cast<std::size_t>(i)],
                                     facet[static_cast<std::size_t>(j)]};
              std::ranges::sort(edge);
              edges.insert(edge);
            }
          }
          for (auto i = 0; i < 4; ++i)
          {
            for (auto j = i + 1; j < 4; ++j)
            {
              for (auto k = j + 1; k < 4; ++k)
              {
                auto triangle = std::array{facet[static_cast<std::size_t>(i)],
                                           facet[static_cast<std::size_t>(j)],
                                           facet[static_cast<std::size_t>(k)]};
                std::ranges::sort(triangle);
                triangles.insert(triangle);
              }
            }
          }
        }
        result[slice] =
            static_cast<Int_precision>(vertices.size()) -
            static_cast<Int_precision>(edges.size()) +
            static_cast<Int_precision>(triangles.size()) -
            static_cast<Int_precision>(facets_by_slice[slice].size());
      }
      return result;
    }

    [[nodiscard]] auto spatial_slices_are_connected() const -> bool
    {
      if (m_simplices.empty()) { return false; }
      auto const facets_by_slice = spacelike_facets_by_slice();
      for (auto const& facets : facets_by_slice)
      {
        if (facets.empty()) { continue; }
        std::map<std::array<VertexId, 3>, std::vector<std::array<VertexId, 4>>>
            triangle_to_facets;
        for (auto const& facet : facets)
        {
          for (auto omitted = 0; omitted < 4; ++omitted)
          {
            std::array<VertexId, 3> triangle{};
            auto                    out = 0;
            for (auto index = 0; index < 4; ++index)
            {
              if (index == omitted) { continue; }
              triangle[static_cast<std::size_t>(out++)] =
                  facet[static_cast<std::size_t>(index)];
            }
            std::ranges::sort(triangle);
            triangle_to_facets[triangle].push_back(facet);
          }
        }

        std::set<std::array<VertexId, 4>>   visited;
        std::queue<std::array<VertexId, 4>> frontier;
        frontier.push(*facets.begin());
        visited.insert(*facets.begin());
        while (!frontier.empty())
        {
          auto const facet = frontier.front();
          frontier.pop();
          for (auto omitted = 0; omitted < 4; ++omitted)
          {
            std::array<VertexId, 3> triangle{};
            auto                    out = 0;
            for (auto index = 0; index < 4; ++index)
            {
              if (index == omitted) { continue; }
              triangle[static_cast<std::size_t>(out++)] =
                  facet[static_cast<std::size_t>(index)];
            }
            std::ranges::sort(triangle);
            for (auto const& neighbor : triangle_to_facets[triangle])
            {
              if (visited.insert(neighbor).second) { frontier.push(neighbor); }
            }
          }
        }
        if (visited.size() != facets.size()) { return false; }
      }
      return true;
    }

    [[nodiscard]] auto simplex_neighbor_graph_connected() const -> bool
    {
      if (m_simplices.empty()) { return true; }
      std::unordered_map<SimplexId, Simplex4D const*> simplex_by_id;
      for (auto const& simplex : m_simplices)
      {
        simplex_by_id.emplace(simplex.id, &simplex);
      }
      std::set<SimplexId>   visited;
      std::queue<SimplexId> frontier;
      frontier.push(m_simplices.front().id);
      visited.insert(m_simplices.front().id);
      while (!frontier.empty())
      {
        auto const current = frontier.front();
        frontier.pop();
        auto const current_it = simplex_by_id.find(current);
        if (current_it == simplex_by_id.end()) { continue; }
        for (auto const& neighbor : current_it->second->neighbors)
        {
          if (!neighbor) { continue; }
          if (visited.insert(*neighbor).second) { frontier.push(*neighbor); }
        }
      }
      return visited.size() == m_simplices.size();
    }

    [[nodiscard]] auto topology_matches_closed_s3_slices() const -> bool
    {
      auto const eulers = derived_slice_euler_characteristics();
      return std::ranges::all_of(eulers,
                                 [](auto const chi) { return chi == 0; }) &&
             (!m_closed_s3_slices || spatial_slices_are_connected());
    }

    void add_vertex(VertexId& next_vertex, Int_precision const time)
    {
      m_vertices.push_back(Vertex4D{next_vertex, time});
      m_vertex_times[next_vertex] = time;
      ++next_vertex;
    }

   public:
    FoliatedTriangulation4() = default;

    explicit FoliatedTriangulation4(Int_precision const timeslices)
    {
      *this = periodic_seed(timeslices);
    }

    FoliatedTriangulation4(Int_precision const timeslices, S4Counts counts,
                           Profile profile)
        : m_timeslices{std::max<Int_precision>(2, timeslices)}
        , m_periodic{true}
        , m_counts{counts}
        , m_proposal_inventory{proposal_inventory_from_counts(m_counts)}
        , m_spatial_profile{std::move(profile)}
        , m_closed_s3_slices{true}
    {
      if (m_spatial_profile.empty())
      {
        m_spatial_profile.assign(static_cast<std::size_t>(m_timeslices), 0);
      }
    }

    [[nodiscard]] static auto from_counts_for_validation(
        Int_precision const timeslices, S4Counts counts, Profile profile)
        -> FoliatedTriangulation4
    {
      return FoliatedTriangulation4{timeslices, counts, std::move(profile)};
    }

    [[nodiscard]] static auto from_checkpoint_state(
        Int_precision const timeslices, S4Counts counts, Profile profile,
        VertexContainer vertices, SimplexContainer simplices,
        bool const three_three_forward) -> FoliatedTriangulation4
    {
      FoliatedTriangulation4 result{timeslices, counts, std::move(profile)};
      result.m_vertices            = std::move(vertices);
      result.m_simplices           = std::move(simplices);
      result.m_three_three_forward = three_three_forward;
      result.rebuild_vertex_time_cache();
      if (!result.m_vertices.empty() || !result.m_simplices.empty())
      {
        result.refresh_derived_state();
      }
      else
      {
        result.m_proposal_inventory =
            proposal_inventory_from_counts(result.m_counts);
      }
      return result;
    }

    [[nodiscard]] static auto periodic_seed(Int_precision const timeslices)
        -> FoliatedTriangulation4
    {
      FoliatedTriangulation4 result;
      result.m_timeslices = std::max<Int_precision>(2, timeslices);
      result.m_periodic   = true;
      result.m_vertices.clear();
      result.m_simplices.clear();
      result.m_vertex_times.clear();

      VertexId                             next_vertex  = 1;
      SimplexId                            next_simplex = 1;
      std::vector<std::array<VertexId, 5>> slice_vertices(
          static_cast<std::size_t>(result.m_timeslices));
      for (auto time = 0; time < result.m_timeslices; ++time)
      {
        for (auto vertex = 0; vertex < 5; ++vertex)
        {
          slice_vertices[static_cast<std::size_t>(time)]
                        [static_cast<std::size_t>(vertex)] = next_vertex;
          result.add_vertex(next_vertex, time);
        }
      }

      for (auto time = 0; time < result.m_timeslices; ++time)
      {
        auto const  next_time = (time + 1) % result.m_timeslices;
        auto const& lower     = slice_vertices[static_cast<std::size_t>(time)];
        auto const& upper = slice_vertices[static_cast<std::size_t>(next_time)];

        auto add_simplex = [&](std::array<VertexId, 5> vertices) {
          Simplex4D simplex;
          simplex.id       = next_simplex++;
          simplex.vertices = vertices;
          simplex.type     = result.classify_simplex(vertices).value();
          result.m_simplices.push_back(simplex);
        };

        if (result.m_timeslices == 2)
        {
          for (auto omitted_base_vertex = 0; omitted_base_vertex < 5;
               ++omitted_base_vertex)
          {
            std::array<int, 4> base{};
            auto               out = 0;
            for (auto vertex = 0; vertex < 5; ++vertex)
            {
              if (vertex != omitted_base_vertex)
              {
                base[static_cast<std::size_t>(out++)] = vertex;
              }
            }

            add_simplex(std::array<VertexId, 5>{
                lower[static_cast<std::size_t>(base[0])],
                lower[static_cast<std::size_t>(base[1])],
                lower[static_cast<std::size_t>(base[2])],
                lower[static_cast<std::size_t>(base[3])],
                upper[static_cast<std::size_t>(base[3])]});
            add_simplex(std::array<VertexId, 5>{
                lower[static_cast<std::size_t>(base[0])],
                lower[static_cast<std::size_t>(base[1])],
                lower[static_cast<std::size_t>(base[2])],
                upper[static_cast<std::size_t>(base[2])],
                upper[static_cast<std::size_t>(base[3])]});
            add_simplex(std::array<VertexId, 5>{
                lower[static_cast<std::size_t>(base[0])],
                lower[static_cast<std::size_t>(base[1])],
                upper[static_cast<std::size_t>(base[1])],
                upper[static_cast<std::size_t>(base[2])],
                upper[static_cast<std::size_t>(base[3])]});
            add_simplex(std::array<VertexId, 5>{
                lower[static_cast<std::size_t>(base[0])],
                upper[static_cast<std::size_t>(base[0])],
                upper[static_cast<std::size_t>(base[1])],
                upper[static_cast<std::size_t>(base[2])],
                upper[static_cast<std::size_t>(base[3])]});
          }
          continue;
        }

        for (auto mask = 1; mask < (1 << 5) - 1; ++mask)
        {
          std::array<VertexId, 5> vertices{};
          auto                    out         = 0;
          auto                    lower_count = 0;
          for (auto vertex = 0; vertex < 5; ++vertex)
          {
            if ((mask & (1 << vertex)) != 0)
            {
              vertices[static_cast<std::size_t>(out++)] =
                  lower[static_cast<std::size_t>(vertex)];
              ++lower_count;
            }
          }
          for (auto vertex = 4; vertex >= 0; --vertex)
          {
            if ((mask & (1 << vertex)) == 0)
            {
              vertices[static_cast<std::size_t>(out++)] =
                  upper[static_cast<std::size_t>(4 - vertex)];
            }
          }
          if (lower_count > 0 && lower_count < 5) { add_simplex(vertices); }
        }
      }

      result.m_closed_s3_slices = true;
      result.refresh_derived_state();
      return result;
    }

    [[nodiscard]] auto timeslices() const -> Int_precision
    {
      return m_timeslices;
    }

    [[nodiscard]] auto periodic() const -> bool { return m_periodic; }

    [[nodiscard]] auto vertices() const -> VertexContainer const&
    {
      return m_vertices;
    }

    [[nodiscard]] auto simplices() const -> SimplexContainer const&
    {
      return m_simplices;
    }

    [[nodiscard]] auto counts() const -> S4Counts { return m_counts; }

    [[nodiscard]] auto proposal_inventory() const -> ProposalInventory4D
    {
      return m_proposal_inventory;
    }

    [[nodiscard]] auto three_three_forward() const -> bool
    {
      return m_three_three_forward;
    }

    [[nodiscard]] auto has_closed_s3_slices() const -> bool
    {
      return !m_simplices.empty() && m_closed_s3_slices &&
             topology_matches_closed_s3_slices();
    }

    [[nodiscard]] auto spatial_topology() const -> std::string_view
    {
      return has_closed_s3_slices() ? "S3" : "non-S3";
    }

    [[nodiscard]] auto spacetime_topology() const -> std::string_view
    {
      if (!has_closed_s3_slices()) { return "unvalidated"; }
      return m_periodic ? "S3xS1" : "S3xI";
    }

    [[nodiscard]] auto slice_euler_characteristics() const
        -> std::vector<Int_precision>
    {
      return derived_slice_euler_characteristics();
    }

    [[nodiscard]] auto spatial_volume_profile() const -> Profile
    {
      return m_spatial_profile;
    }

    [[nodiscard]] auto centered_spatial_volume_profile() const -> Profile
    {
      auto profile = m_spatial_profile;
      if (profile.empty()) { return profile; }
      auto const peak = static_cast<std::size_t>(
          std::distance(profile.begin(), std::ranges::max_element(profile)));
      auto const center = profile.size() / 2;
      std::rotate(profile.begin(),
                  profile.begin() +
                      static_cast<std::ptrdiff_t>(
                          (peak + profile.size() - center) % profile.size()),
                  profile.end());
      return profile;
    }

    [[nodiscard]] auto max_vertex_order() const -> Int_precision
    {
      std::map<VertexId, Int_precision> orders;
      for (auto const& simplex : m_simplices)
      {
        for (auto const vertex : simplex.vertices) { ++orders[vertex]; }
      }
      auto const max_order = std::ranges::max_element(
          orders, {}, [](auto const& pair) { return pair.second; });
      auto const move_growth = std::max<Int_precision>(
          0, m_counts.N4 - static_cast<Int_precision>(m_simplices.size()));
      return max_order == orders.end() ? move_growth
                                       : max_order->second + move_growth;
    }

    [[nodiscard]] auto vertex_order_distribution() const
        -> std::map<Int_precision, Int_precision>
    {
      std::map<VertexId, Int_precision> orders;
      for (auto const& simplex : m_simplices)
      {
        for (auto const vertex : simplex.vertices) { ++orders[vertex]; }
      }
      std::map<Int_precision, Int_precision> distribution;
      for (auto const& [_, order] : orders) { ++distribution[order]; }
      return distribution;
    }

    [[nodiscard]] auto occupied_temporal_width() const -> Int_precision
    {
      return static_cast<Int_precision>(std::ranges::count_if(
          m_spatial_profile, [](auto const volume) { return volume > 0; }));
    }

    [[nodiscard]] auto slice_to_slice_roughness() const -> long double
    {
      if (m_spatial_profile.size() < 2) { return 0.0L; }
      long double roughness = 0.0L;
      for (std::size_t index = 0; index < m_spatial_profile.size(); ++index)
      {
        auto const next = (index + 1) % m_spatial_profile.size();
        roughness += std::abs(static_cast<long double>(
            m_spatial_profile[index] - m_spatial_profile[next]));
      }
      return roughness;
    }

    [[nodiscard]] auto inverse_participation_ratio() const -> long double
    {
      auto const total = std::accumulate(m_spatial_profile.begin(),
                                         m_spatial_profile.end(), 0.0L);
      if (total == 0.0L) { return 0.0L; }
      auto square_sum = 0.0L;
      for (auto const volume : m_spatial_profile)
      {
        square_sum +=
            static_cast<long double>(volume) * static_cast<long double>(volume);
      }
      return square_sum / (total * total);
    }

    [[nodiscard]] auto alternating_slice_order_parameter() const -> long double
    {
      auto const total = std::accumulate(m_spatial_profile.begin(),
                                         m_spatial_profile.end(), 0.0L);
      if (total == 0.0L) { return 0.0L; }
      auto alternating = 0.0L;
      for (std::size_t index = 0; index < m_spatial_profile.size(); ++index)
      {
        alternating += (index % 2 == 0 ? 1.0L : -1.0L) *
                       static_cast<long double>(m_spatial_profile[index]);
      }
      return alternating / total;
    }

    [[nodiscard]] static auto move_count_delta(move_tracker::MoveType4D move)
        -> S4Counts
    {
      return move_descriptor_4d(move).delta;
    }

    [[nodiscard]] auto candidate_multiplicity(
        move_tracker::MoveType4D const move) const -> Int_precision
    {
      using move_tracker::MoveType4D;
      switch (move)
      {
        case MoveType4D::TWO_FOUR:
          return static_cast<Int_precision>(enumerate_two_four_sites().size());
        case MoveType4D::FOUR_TWO:
          return static_cast<Int_precision>(enumerate_four_two_sites().size());
        case MoveType4D::THREE_THREE:
          return static_cast<Int_precision>(
              enumerate_three_three_sites().size());
        case MoveType4D::FOUR_SIX:
          return static_cast<Int_precision>(enumerate_four_six_sites().size());
        case MoveType4D::SIX_FOUR:
          return static_cast<Int_precision>(enumerate_six_four_sites().size());
        case MoveType4D::TWO_EIGHT:
          return static_cast<Int_precision>(enumerate_two_eight_sites().size());
        case MoveType4D::EIGHT_TWO:
          return static_cast<Int_precision>(
              enumerate_eight_two_sites().size());
        default: return 0;
      }
    }

    [[nodiscard]] auto is_applicable(move_tracker::MoveType4D const move) const
        -> bool
    {
      return candidate_multiplicity(move) > 0;
    }

    [[nodiscard]] auto apply_move(move_tracker::MoveType4D const move,
                                  std::size_t const site_index = 0) -> bool
    {
      using move_tracker::MoveType4D;
      if (move == MoveType4D::TWO_FOUR)
      {
        auto const sites = enumerate_two_four_sites();
        if (site_index >= sites.size()) { return false; }
        auto const&                          site = sites[site_index];
        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(4);
        for (auto omitted_vertex : site.shared_facet)
        {
          std::array<VertexId, 5> vertices{};
          auto                    out = 0;
          for (auto const vertex : site.shared_facet)
          {
            if (vertex != omitted_vertex)
            {
              vertices[static_cast<std::size_t>(out++)] = vertex;
            }
          }
          vertices[static_cast<std::size_t>(out++)] = site.new_edge[0];
          vertices[static_cast<std::size_t>(out++)] = site.new_edge[1];
          replacement.push_back(vertices);
        }
        return replace_simplices(
            std::set<SimplexId>{site.simplex_ids[0], site.simplex_ids[1]},
            replacement);
      }
      if (move == MoveType4D::FOUR_TWO)
      {
        auto const sites = enumerate_four_two_sites();
        if (site_index >= sites.size()) { return false; }
        auto const& site              = sites[site_index];
        auto        first_replacement = std::array<VertexId, 5>{
            site.replacement_facet[0], site.replacement_facet[1],
            site.replacement_facet[2], site.replacement_facet[3], site.edge[0]};
        auto second_replacement = std::array<VertexId, 5>{
            site.replacement_facet[0], site.replacement_facet[1],
            site.replacement_facet[2], site.replacement_facet[3], site.edge[1]};
        auto const removed_ids = std::set<SimplexId>(site.simplex_ids.begin(),
                                                     site.simplex_ids.end());
        return replace_simplices(
            removed_ids, std::vector{first_replacement, second_replacement});
      }
      if (move == MoveType4D::THREE_THREE)
      {
        auto const sites = enumerate_three_three_sites();
        if (site_index >= sites.size()) { return false; }
        auto const& site = sites[site_index];
        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(3);
        for (auto omitted_vertex : site.triangle)
        {
          std::array<VertexId, 5> vertices{};
          auto                    out = 0;
          for (auto const vertex : site.replacement_triangle)
          {
            vertices[static_cast<std::size_t>(out++)] = vertex;
          }
          for (auto const vertex : site.triangle)
          {
            if (vertex != omitted_vertex)
            {
              vertices[static_cast<std::size_t>(out++)] = vertex;
            }
          }
          replacement.push_back(vertices);
        }
        auto const removed_ids = std::set<SimplexId>(
            site.simplex_ids.begin(), site.simplex_ids.end());
        return replace_simplices(removed_ids, replacement);
      }
      if (move == MoveType4D::FOUR_SIX)
      {
        auto const sites = enumerate_four_six_sites();
        if (site_index >= sites.size()) { return false; }
        auto const& site = sites[site_index];
        std::set<VertexId> apices;
        for (auto const id : site.simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr) { return false; }
          for (auto const vertex : simplex->vertices)
          {
            if (!contains_vertex(site.triangle, vertex) &&
                !contains_vertex(site.new_edge, vertex))
            {
              apices.insert(vertex);
            }
          }
        }
        if (apices.size() != 2) { return false; }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(6);
        for (auto const apex : apices)
        {
          for (auto const omitted_vertex : site.triangle)
          {
            std::array<VertexId, 5> vertices{};
            auto                    out = 0;
            vertices[static_cast<std::size_t>(out++)] = site.new_edge[0];
            vertices[static_cast<std::size_t>(out++)] = site.new_edge[1];
            vertices[static_cast<std::size_t>(out++)] = apex;
            for (auto const vertex : site.triangle)
            {
              if (vertex != omitted_vertex)
              {
                vertices[static_cast<std::size_t>(out++)] = vertex;
              }
            }
            replacement.push_back(vertices);
          }
        }
        auto const removed_ids = std::set<SimplexId>(
            site.simplex_ids.begin(), site.simplex_ids.end());
        return replace_simplices(removed_ids, replacement);
      }
      if (move == MoveType4D::SIX_FOUR)
      {
        auto const sites = enumerate_six_four_sites();
        if (site_index >= sites.size()) { return false; }
        auto const& site = sites[site_index];
        std::set<VertexId> apices;
        for (auto const id : site.simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr) { return false; }
          for (auto const vertex : simplex->vertices)
          {
            if (!contains_vertex(site.replacement_triangle, vertex) &&
                !contains_vertex(site.edge, vertex))
            {
              apices.insert(vertex);
            }
          }
        }
        if (apices.size() != 2) { return false; }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(4);
        for (auto const endpoint : site.edge)
        {
          for (auto const apex : apices)
          {
            replacement.push_back(std::array<VertexId, 5>{
                site.replacement_triangle[0], site.replacement_triangle[1],
                site.replacement_triangle[2], endpoint, apex});
          }
        }
        auto const removed_ids = std::set<SimplexId>(
            site.simplex_ids.begin(), site.simplex_ids.end());
        return replace_simplices(removed_ids, replacement);
      }
      if (move == MoveType4D::TWO_EIGHT)
      {
        auto const sites = enumerate_two_eight_sites();
        if (site_index >= sites.size()) { return false; }
        auto const& site       = sites[site_index];
        auto const  new_vertex = next_vertex_id();
        std::set<VertexId> apices;
        for (auto const id : site.simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr) { return false; }
          for (auto const vertex : simplex->vertices)
          {
            if (!contains_vertex(site.shared_facet, vertex))
            {
              apices.insert(vertex);
            }
          }
        }
        if (apices.size() != 2) { return false; }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(8);
        for (auto const apex : apices)
        {
          for (auto const omitted_vertex : site.shared_facet)
          {
            std::array<VertexId, 5> vertices{};
            auto                    out = 0;
            vertices[static_cast<std::size_t>(out++)] = new_vertex;
            vertices[static_cast<std::size_t>(out++)] = apex;
            for (auto const vertex : site.shared_facet)
            {
              if (vertex != omitted_vertex)
              {
                vertices[static_cast<std::size_t>(out++)] = vertex;
              }
            }
            replacement.push_back(vertices);
          }
        }
        auto const removed_ids =
            std::set<SimplexId>{site.simplex_ids[0], site.simplex_ids[1]};
        return replace_simplices_and_vertices(
            removed_ids, replacement,
            std::vector<Vertex4D>{Vertex4D{new_vertex, site.vertex_time}});
      }
      if (move == MoveType4D::EIGHT_TWO)
      {
        auto const sites = enumerate_eight_two_sites();
        if (site_index >= sites.size()) { return false; }
        auto const& site = sites[site_index];
        std::set<VertexId> apices;
        for (auto const id : site.simplex_ids)
        {
          auto const* simplex = simplex_by_id(id);
          if (simplex == nullptr) { return false; }
          for (auto const vertex : simplex->vertices)
          {
            if (!contains_vertex(site.replacement_facet, vertex) &&
                vertex != site.vertex)
            {
              apices.insert(vertex);
            }
          }
        }
        if (apices.size() != 2) { return false; }

        std::vector<std::array<VertexId, 5>> replacement;
        replacement.reserve(2);
        for (auto const apex : apices)
        {
          replacement.push_back(std::array<VertexId, 5>{
              site.replacement_facet[0], site.replacement_facet[1],
              site.replacement_facet[2], site.replacement_facet[3], apex});
        }
        auto const removed_ids = std::set<SimplexId>(
            site.simplex_ids.begin(), site.simplex_ids.end());
        return replace_simplices_and_vertices(
            removed_ids, replacement, std::vector<Vertex4D>{},
            std::set<VertexId>{site.vertex});
      }
      return false;
    }

    [[nodiscard]] auto validate() const -> ValidationReport
    {
      ValidationReport report;
      report.standard_cdt_candidate = true;
      if (m_timeslices < 2)
      {
        report.errors.emplace_back("At least two timeslices are required.");
      }
      if (m_counts.N4 !=
          m_counts.N41 + m_counts.N32 + m_counts.N23 + m_counts.N14)
      {
        report.errors.emplace_back("N4 does not match the simplex type sum.");
      }
      if (m_counts.N0 < 0 || m_counts.N1 < 0 || m_counts.N2 < 0 ||
          m_counts.N3 < 0 || m_counts.N4 < 0 || m_counts.N41 < 0 ||
          m_counts.N32 < 0 || m_counts.N23 < 0 || m_counts.N14 < 0)
      {
        report.errors.emplace_back("Negative simplex count found.");
      }
      if (m_vertices.empty() || m_simplices.empty())
      {
        report.standard_cdt_candidate = false;
        report.errors.emplace_back(
            "A standard 4D CDT candidate requires an explicit simplex complex.");
      }
      if (!m_periodic)
      {
        report.standard_cdt_candidate = false;
        report.errors.emplace_back(
            "Standard CDT candidate requires periodic time.");
      }
      if (!m_closed_s3_slices)
      {
        report.standard_cdt_candidate = false;
        report.errors.emplace_back(
            "Spatial slices are not marked as closed S3.");
      }
      if (!topology_matches_closed_s3_slices())
      {
        report.standard_cdt_candidate = false;
        report.errors.emplace_back(
            "Spatial slices are not validated as connected S3 slices.");
      }
      for (auto const chi : slice_euler_characteristics())
      {
        if (chi != 0)
        {
          report.standard_cdt_candidate = false;
          report.errors.emplace_back(
              "A spatial slice does not have S3 Euler characteristic.");
          break;
        }
      }
      if (m_spatial_profile.size() != static_cast<std::size_t>(m_timeslices))
      {
        report.errors.emplace_back(
            "Spatial profile does not match timeslice count.");
      }
      if (std::ranges::any_of(m_spatial_profile,
                              [](auto const volume) { return volume < 0; }))
      {
        report.errors.emplace_back("Spatial profile contains negative volume.");
      }
      if (m_vertex_times.size() != m_vertices.size())
      {
        report.errors.emplace_back(
            "Vertex-time cache does not match vertices.");
      }
      for (auto const& vertex : m_vertices)
      {
        auto const cached = m_vertex_times.find(vertex.id);
        if (cached == m_vertex_times.end() || cached->second != vertex.time)
        {
          report.errors.emplace_back("Vertex-time cache is stale.");
          break;
        }
      }
      if (m_proposal_inventory.two_four_sites < 0 ||
          m_proposal_inventory.four_two_sites < 0 ||
          m_proposal_inventory.three_three_sites < 0 ||
          m_proposal_inventory.four_six_sites < 0 ||
          m_proposal_inventory.six_four_sites < 0 ||
          m_proposal_inventory.two_eight_sites < 0 ||
          m_proposal_inventory.eight_two_sites < 0)
      {
        report.errors.emplace_back("Negative proposal multiplicity found.");
      }

      std::set<std::array<VertexId, 5>>               simplex_keys;
      std::map<std::array<VertexId, 4>, int>          facet_incidence;
      std::unordered_map<SimplexId, Simplex4D const*> simplex_by_id;
      for (auto const& simplex : m_simplices)
      {
        simplex_by_id.emplace(simplex.id, &simplex);
      }

      for (auto const& simplex : m_simplices)
      {
        auto vertices = sorted_vertices(simplex.vertices);
        if (!simplex_keys.insert(vertices).second)
        {
          report.errors.emplace_back("Duplicate 4-simplex found.");
        }
        if (std::set<VertexId>(simplex.vertices.begin(), simplex.vertices.end())
                .size() != 5)
        {
          report.errors.emplace_back("A 4-simplex has duplicate vertices.");
        }
        if (std::ranges::any_of(simplex.vertices, [&](auto const vertex) {
              return !m_vertex_times.contains(vertex);
            }))
        {
          report.errors.emplace_back(
              "A 4-simplex references a missing vertex.");
          continue;
        }
        auto const expected_type = classify_simplex(simplex.vertices);
        if (!expected_type || *expected_type != simplex.type)
        {
          report.errors.emplace_back("A 4-simplex has invalid causal type.");
        }
        for (auto index = 0; index < 5; ++index)
        {
          auto const facet = facet_vertices(simplex, index);
          ++facet_incidence[facet];
          auto const neighbor = simplex.neighbors[static_cast<size_t>(index)];
          if (!neighbor) { continue; }
          auto const neighbor_it = simplex_by_id.find(*neighbor);
          if (neighbor_it == simplex_by_id.end())
          {
            report.errors.emplace_back("Neighbor simplex ID does not exist.");
            continue;
          }
          auto const& neighbor_simplex = *neighbor_it->second;
          auto const  reciprocal       = std::ranges::any_of(
              neighbor_simplex.neighbors, [&](auto const& maybe_neighbor) {
                return maybe_neighbor && *maybe_neighbor == simplex.id;
              });
          if (!reciprocal)
          {
            report.errors.emplace_back(
                "Neighbor relationship is not reciprocal.");
          }
        }
      }

      if (m_closed_s3_slices && !simplex_neighbor_graph_connected())
      {
        report.errors.emplace_back("Simplex neighbor graph is disconnected.");
      }

      if (m_periodic)
      {
        for (auto const& [_, incidence] : facet_incidence)
        {
          if (incidence != 2)
          {
            report.errors.emplace_back(
                "Periodic triangulation has an unintended boundary.");
            break;
          }
        }
      }

      return report;
    }

    [[nodiscard]] auto is_valid() const -> bool { return validate().valid(); }

    [[nodiscard]] auto time_reversed() const -> FoliatedTriangulation4
    {
      auto reversed = *this;
      for (auto& vertex : reversed.m_vertices)
      {
        vertex.time = (m_timeslices - vertex.time) % m_timeslices;
      }
      reversed.rebuild_vertex_time_cache();
      for (auto& simplex : reversed.m_simplices)
      {
        switch (simplex.type)
        {
          case SimplexType4D::FOUR_ONE:
            simplex.type = SimplexType4D::ONE_FOUR;
            break;
          case SimplexType4D::ONE_FOUR:
            simplex.type = SimplexType4D::FOUR_ONE;
            break;
          case SimplexType4D::THREE_TWO:
            simplex.type = SimplexType4D::TWO_THREE;
            break;
          case SimplexType4D::TWO_THREE:
            simplex.type = SimplexType4D::THREE_TWO;
            break;
        }
      }
      reversed.m_three_three_forward = !m_three_three_forward;
      if (!reversed.m_simplices.empty())
      {
        reversed.refresh_derived_state();
        return reversed;
      }
      std::swap(reversed.m_counts.N41, reversed.m_counts.N14);
      std::swap(reversed.m_counts.N32, reversed.m_counts.N23);
      reversed.m_proposal_inventory = ProposalInventory4D{};
      if (m_spatial_profile.size() == static_cast<std::size_t>(m_timeslices))
      {
        Profile    mapped(m_spatial_profile.size(), 0);
        auto const slices = static_cast<std::size_t>(m_timeslices);
        for (std::size_t index = 0; index < slices; ++index)
        {
          mapped[(slices - index) % slices] = m_spatial_profile[index];
        }
        reversed.m_spatial_profile = std::move(mapped);
      }
      return reversed;
    }

    [[nodiscard]] auto canonical_hash() const -> std::string
    {
      std::ostringstream stream;
      stream << "T=" << m_timeslices << ";P=" << m_periodic
             << ";F=" << m_three_three_forward << ";";
      stream << m_counts.N0 << ',' << m_counts.N1 << ',' << m_counts.N2 << ','
             << m_counts.N3 << ',' << m_counts.N4 << ',' << m_counts.N41 << ','
             << m_counts.N32 << ',' << m_counts.N23 << ',' << m_counts.N14
             << ";C=" << m_counts.class_resolved.has_value();
      if (m_counts.class_resolved)
      {
        auto const& class_counts = *m_counts.class_resolved;
        stream << ',' << class_counts.spatial_tetrahedra << ','
               << class_counts.timelike_edges << ','
               << class_counts.mixed_triangles << ','
               << class_counts.timelike_tetrahedra;
      }
      else { stream << ",0,0,0,0"; }
      stream << ";V=";
      for (auto const volume : m_spatial_profile) { stream << volume << ','; }
      stream << ";VT=";
      std::vector<std::pair<VertexId, Int_precision>> vertices;
      vertices.reserve(m_vertices.size());
      for (auto const& vertex : m_vertices)
      {
        vertices.push_back({vertex.id, vertex.time});
      }
      std::ranges::sort(vertices);
      for (auto const& [id, time] : vertices)
      {
        stream << id << '@' << time << ',';
      }
      stream << ";S=";
      std::vector<std::array<VertexId, 5>> simplices;
      simplices.reserve(m_simplices.size());
      for (auto const& simplex : m_simplices)
      {
        simplices.push_back(sorted_vertices(simplex.vertices));
      }
      std::ranges::sort(simplices);
      for (auto const& simplex : simplices)
      {
        for (auto const vertex : simplex) { stream << vertex << '.'; }
        stream << static_cast<int>(classify_simplex(simplex).value_or(
                      SimplexType4D::FOUR_ONE))
               << ',';
      }
      return stream.str();
    }
  };
}  // namespace cdt::four_d

#endif  // CDT_PLUSPLUS_FOLIATED_TRIANGULATION_4_HPP
