// Boost.Geometry
// Unit Test Helper

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_TEST_BUFFER_CSV_HPP
#define BOOST_GEOMETRY_TEST_BUFFER_CSV_HPP

#include <fstream>
#include <iomanip>

#include "debug_buffer_info.hpp"
#include <boost/geometry/algorithms/detail/overlay/debug_turn_info.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>


namespace detail
{

class buffer_visitor_csv
{
public:
    buffer_visitor_csv(const std::string &path)
        : m_path(path)
    {
    }

    template <typename PieceCollection>
    inline void apply(PieceCollection const& collection, int phase)
    {
        if (phase == 0)
        {
            write_pieces(collection.m_pieces, collection.offsetted_rings);
        }
        if (phase == 1)
        {
            write_turns(collection.m_turns);
        }
    }

    template <typename Geometry, typename GeometryBuffer>
    inline void write_input_output(Geometry const& input, GeometryBuffer const& buffer)
    {
        std::ofstream stream(m_path + "_geometries.csv");
        if (! stream.good())
        {
            return;
        }

        stream << std::setprecision(16)
            << "id;role;wkt" << std::endl
            << "1;input;" << bg::wkt(input) << std::endl
            << "2;output;" << bg::wkt(buffer) << std::endl;
    }

private:

    template <typename Turns>
    inline void write_turns(Turns const& turns)
    {
        std::ofstream stream(m_path + "_turns.csv");
        if (! stream.good())
        {
            return;
        }
        std::size_t id = 1;
        stream << std::setprecision(16) << "id;wkt;piece_indices;method;operations;cluster_id;traversable;blocked" << std::endl;
        for (const auto &turn : turns)
        {
            stream << id++ << ";"
                   << bg::wkt(turn.point) << ";"
                   << turn.operations[0].piece_index << "/" << turn.operations[1].piece_index << ";"
                   << bg::method_char(turn.method) << ";"
                   << bg::operation_char(turn.operations[0].operation) << "/" << bg::operation_char(turn.operations[1].operation) << ";"
                   << turn.cluster_id << ";"
                   << turn.is_turn_traversable << ";"
                   << turn.blocked() << ";"
                   << std::endl;
        }
    }

    template <typename Pieces, typename OffsettedRings>
    inline void write_pieces(Pieces const& pieces, OffsettedRings const& offsetted_rings)
    {
        std::ofstream stream(m_path + "_pieces.csv");
        if (! stream.good())
        {
            return;
        }

        std::size_t id = 1;
        stream << std::setprecision(16) << "id;wkt;index;start;end;type;segment_indices" << std::endl;
        for (auto const& piece : pieces)
        {
            bg::segment_identifier const& seg_id = piece.first_seg_id;
            if (seg_id.segment_index < 0)
            {
                continue;
            }
            // NOTE: ring is returned by copy here
            auto const corner = piece.m_piece_border.get_full_ring();
            if (corner.empty())
            {
                continue;
            }

            auto const& ring = offsetted_rings[seg_id.multi_index];

            stream << id++ << ";"
                   << bg::wkt(corner) << ";"
                   << piece.index << ";"
                   << piece.is_flat_start << ";"
                   << piece.is_flat_end << ";"
                   << bg::debug::piece_type_char(piece.type) << ";"
                   << piece.first_seg_id.segment_index << ".." << piece.beyond_last_segment_index - 1
                   << std::endl;
        }
    }

    std::string m_path;
};

}

#endif
