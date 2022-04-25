// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test Helper

// Copyright (c) 2010-2015 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2020-2021.
// Modifications copyright (c) 2020-2021 Oracle and/or its affiliates.
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_GEOMETRY_TEST_BUFFER_SVG_HPP
#define BOOST_GEOMETRY_TEST_BUFFER_SVG_HPP

#include <fstream>
#include <sstream>

// Uncomment next lines if you want to have a zoomed view
//#define BOOST_GEOMETRY_BUFFER_TEST_SVG_USE_ALTERNATE_BOX

// If possible define box before including this unit with the right view
#ifdef BOOST_GEOMETRY_BUFFER_TEST_SVG_USE_ALTERNATE_BOX
#  ifndef BOOST_GEOMETRY_BUFFER_TEST_SVG_ALTERNATE_BOX
#    define BOOST_GEOMETRY_BUFFER_TEST_SVG_ALTERNATE_BOX "BOX(0 0,100 100)"
#  endif
#endif

#include <boost/geometry/io/svg/svg_mapper.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/intersection.hpp>

namespace detail
{

template <typename Ring, typename cstag = typename bg::cs_tag<Ring>::type>
struct get_labelpoint
{
    using point_type = typename bg::point_type<Ring>::type;

    template <typename Piece>
    static point_type apply(Ring const& , Piece const& piece)
    {
        return piece.m_label_point;
    }
};

template <typename Ring>
struct get_labelpoint<Ring, bg::cartesian_tag>
{
    using point_type = typename bg::point_type<Ring>::type;

    template <typename Piece>
    static point_type apply(Ring const& ring, Piece const& piece)
    {
        // Centroid is currently only available for cartesian
        return ring.empty()
                ? piece.m_label_point
                : bg::return_centroid<point_type>(ring);
    }
};

}

template <typename Ring, typename Piece>
inline typename bg::point_type<Ring>::type
get_labelpoint(Ring const& ring, Piece const& piece)
{
    if ((piece.type == bg::strategy::buffer::buffered_concave
         || piece.type == bg::strategy::buffer::buffered_flat_end)
        && ring.size() >= 2u)
    {
        // Return a point between the first two points on the ring
        typename bg::point_type<Ring>::type result;
        bg::set<0>(result, (bg::get<0>(ring[0]) + bg::get<0>(ring[1])) / 2.0);
        bg::set<1>(result, (bg::get<1>(ring[0]) + bg::get<1>(ring[1])) / 2.0);
        return result;
    }
    else
    {
        // Return the piece's labelpoint or the centroid
        return detail::get_labelpoint<Ring>::apply(ring, piece);
    }
}

inline char piece_type_char(bg::strategy::buffer::piece_type const& type)
{
    using namespace bg::strategy::buffer;
    switch(type)
    {
        case buffered_segment : return 's';
        case buffered_join : return 'j';
        case buffered_round_end : return 'r';
        case buffered_flat_end : return 'f';
        case buffered_point : return 'p';
        case buffered_concave : return 'c';
        default : return '?';
    }
}

template <typename SvgMapper, typename Box>
class svg_visitor
{
public :
    svg_visitor(SvgMapper& mapper)
        : m_mapper(mapper)
        , m_zoom(false)
    {
        bg::assign_inverse(m_alternate_box);
    }

    void set_alternate_box(Box const& box)
    {
        m_alternate_box = box;
        m_zoom = true;
    }

    template <typename PieceCollection>
    inline void apply(PieceCollection const& collection, int phase)
    {
        // Comment next return if you want to see pieces, turns, etc.
        return;

        if(phase == 0)
        {
            map_pieces(collection.m_pieces, collection.offsetted_rings, true, true);
        }
        if (phase == 1)
        {
            map_turns(collection.m_turns, true, false);
        }
        if (phase == 2 && ! m_zoom)
        {
//        map_traversed_rings(collection.traversed_rings);
//        map_offsetted_rings(collection.offsetted_rings);
        }
    }

private :
    class si
    {
    private :
        bg::segment_identifier m_id;

    public :
        inline si(bg::segment_identifier const& id)
            : m_id(id)
        {}

        template <typename Char, typename Traits>
        inline friend std::basic_ostream<Char, Traits>& operator<<(
                std::basic_ostream<Char, Traits>& os,
                si const& s)
        {
            os << s.m_id.multi_index << "." << s.m_id.segment_index;
            return os;
        }
    };

    template <typename Turns>
    inline void map_turns(Turns const& turns, bool label_good_turns, bool label_wrong_turns)
    {
        namespace bgdb = boost::geometry::detail::buffer;
        using turn_type = typename boost::range_value<Turns const>::type;
        using point_type = typename turn_type::point_type;

        std::map<point_type, int, bg::less<point_type> > offsets;

        for (auto it = boost::begin(turns); it != boost::end(turns); ++it)
        {
            if (m_zoom && bg::disjoint(it->point, m_alternate_box))
            {
                continue;
            }

            bool is_good = true;
            std::string fill = "fill:rgb(0,255,0);";
            if (! it->is_turn_traversable)
            {
                fill = "fill:rgb(255,0,0);";
                is_good = false;
            }
            if (it->blocked())
            {
                fill = "fill:rgb(128,128,128);";
                is_good = false;
            }

            fill += "fill-opacity:0.7;";

            m_mapper.map(it->point, fill, 4);

            if ((label_good_turns && is_good) || (label_wrong_turns && ! is_good))
            {
                std::ostringstream out;
                out << it->turn_index;
                if (it->cluster_id >= 0)
                {
                   out << " ("  << it->cluster_id << ")";
                }
                out
                    << " " << it->operations[0].piece_index << "/" << it->operations[1].piece_index
                    << " " << si(it->operations[0].seg_id) << "/" << si(it->operations[1].seg_id)

    //              If you want to see travel information
                    << std::endl
                    << " nxt " << it->operations[0].enriched.get_next_turn_index()
                    << "/" << it->operations[1].enriched.get_next_turn_index()
                    //<< " frac " << it->operations[0].fraction

    //                If you want to see point coordinates (e.g. to find duplicates)
                    << std::endl << std::setprecision(16) << bg::wkt(it->point)

                    << std::endl;
                out << " " << bg::method_char(it->method)
                    << ":" << bg::operation_char(it->operations[0].operation)
                    << "/" << bg::operation_char(it->operations[1].operation);
                out << " "
                    << (it->is_turn_traversable ? "" : "w")
                    ;

                offsets[it->point] += 10;
                int offset = offsets[it->point];

                m_mapper.text(it->point, out.str(), "fill:rgb(0,0,0);font-family='Arial';font-size:9px;", 5, offset);

                offsets[it->point] += 25;
            }
        }
    }

    template <typename Pieces, typename OffsettedRings>
    inline void map_pieces(Pieces const& pieces,
                OffsettedRings const& offsetted_rings,
                bool do_pieces, bool do_indices)
    {
        using piece_type = typename boost::range_value<Pieces const>::type ;
        using ring_type = typename boost::range_value<OffsettedRings const>::type;

        for (auto it = boost::begin(pieces); it != boost::end(pieces); ++it)
        {
            const piece_type& piece = *it;
            bg::segment_identifier seg_id = piece.first_seg_id;
            if (seg_id.segment_index < 0)
            {
                continue;
            }

            ring_type const& ring = offsetted_rings[seg_id.multi_index];

#if 0 // Does not compile (SVG is not enabled by default)
            if (m_zoom && bg::disjoint(corner, m_alternate_box))
            {
                continue;
            }
#endif

            // NOTE: ring is returned by copy here
            auto const corner = piece.m_piece_border.get_full_ring();

            if (m_zoom && do_pieces)
            {
                try
                {
                    std::string style = "opacity:0.3;stroke:rgb(0,0,0);stroke-width:1;";
                    bg::model::multi_polygon
                        <
                            bg::model::polygon<typename bg::point_type<Box>::type>
                        > clipped;
                    bg::intersection(ring, m_alternate_box, clipped);
                    m_mapper.map(clipped,
                        piece.type == bg::strategy::buffer::buffered_segment
                        ? style + "fill:rgb(255,128,0);"
                        : style + "fill:rgb(255,0,0);");
                }
                catch (...)
                {
                    std::cerr << "Error for piece " << piece.index << std::endl;
                }
            }
            else if (do_pieces && ! corner.empty())
            {
                std::string style = "opacity:0.3;stroke:rgb(0,0,0);stroke-width:1;";
                m_mapper.map(corner,
                    piece.type == bg::strategy::buffer::buffered_segment
                    ? style + "fill:rgb(255,128,0);"
                    : style + "fill:rgb(255,0,0);");
            }

            if (do_indices)
            {
                // Label starting piece_index / segment_index

                std::ostringstream out;
                out << piece.index
                    << (piece.is_flat_start ? " FS" : "")
                    << (piece.is_flat_end ? " FE" : "")
                    << " (" << piece_type_char(piece.type) << ") "
                    << piece.first_seg_id.segment_index
                    << ".." << piece.beyond_last_segment_index - 1
                       ;

                m_mapper.text(get_labelpoint(corner, piece), out.str(),
                    "fill:rgb(255,0,0);font-family='Arial';font-size:10px;", 5, 5);
            }
        }
    }

    template <typename TraversedRings>
    inline void map_traversed_rings(TraversedRings const& traversed_rings)
    {
        for (auto it = boost::begin(traversed_rings); it != boost::end(traversed_rings); ++it)
        {
            m_mapper.map(*it, "opacity:0.4;fill:none;stroke:rgb(0,255,0);stroke-width:2");
        }
    }

    template <typename OffsettedRings>
    inline void map_offsetted_rings(OffsettedRings const& offsetted_rings)
    {
        for (auto it = boost::begin(offsetted_rings); it != boost::end(offsetted_rings); ++it)
        {
            if (it->discarded())
            {
                m_mapper.map(*it, "opacity:0.4;fill:none;stroke:rgb(255,0,0);stroke-width:2");
            }
            else
            {
                m_mapper.map(*it, "opacity:0.4;fill:none;stroke:rgb(0,0,255);stroke-width:2");
            }
        }
    }


    SvgMapper& m_mapper;
    Box m_alternate_box;
    bool m_zoom;

};

template <typename Point>
class buffer_svg_mapper
{
public :

    buffer_svg_mapper(std::string const& casename)
        : m_casename(casename)
    {
        bg::assign_inverse(m_alternate_box);
    }

    template <typename Mapper, typename Visitor, typename Envelope, typename DistanceType>
    void prepare(Mapper& mapper, Visitor& visitor, Envelope const& envelope,
                 const DistanceType& box_buffer_distance)
    {
#ifdef BOOST_GEOMETRY_BUFFER_TEST_SVG_USE_ALTERNATE_BOX
        // Create a zoomed-in view
        bg::model::box<Point> alternate_box;
        bg::read_wkt(BOOST_GEOMETRY_BUFFER_TEST_SVG_ALTERNATE_BOX, alternate_box);
        mapper.add(alternate_box);

        // Take care non-visible elements are skipped
        visitor.set_alternate_box(alternate_box);
        set_alternate_box(alternate_box);
#else
        bg::model::box<Point> box = envelope;
        bg::buffer(box, box, box_buffer_distance);
        mapper.add(box);
#endif

        boost::ignore_unused(visitor);
    }

    void set_alternate_box(bg::model::box<Point> const& box)
    {
        m_alternate_box = box;
        m_zoom = true;
    }

    template <typename Mapper, typename Geometry, typename GeometryBuffer>
    void map_input_output(Mapper& mapper, Geometry const& geometry,
            GeometryBuffer const& buffered, bool negative)
    {
        bool const areal = bg::util::is_areal<Geometry>::value;

        if (m_zoom)
        {
            map_io_zoomed(mapper, geometry, buffered, negative, areal);
        }
        else
        {
            map_io(mapper, geometry, buffered, negative, areal);
        }
    }

    template <typename Mapper, typename Geometry, typename Strategy, typename RescalePolicy>
    void map_self_ips(Mapper& mapper, Geometry const& geometry, Strategy const& strategy, RescalePolicy const& rescale_policy)
    {
        using turn_info = bg::detail::overlay::turn_info
        <
            Point,
            typename bg::detail::segment_ratio_type<Point, RescalePolicy>::type
        >;

        std::vector<turn_info> turns;

        bg::detail::self_get_turn_points::no_interrupt_policy policy;
        bg::self_turns
            <
                bg::detail::overlay::assign_null_policy
            >(geometry, strategy, rescale_policy, turns, policy);

        for (turn_info const& turn : turns)
        {
            mapper.map(turn.point, "fill:rgb(255,128,0);stroke:rgb(0,0,100);stroke-width:1", 3);
        }
    }

private :

    template <typename Mapper, typename Geometry, typename GeometryBuffer>
    void map_io(Mapper& mapper, Geometry const& geometry,
            GeometryBuffer const& buffered, bool negative, bool areal)
    {
        // Map input geometry in green
        if (areal)
        {
            mapper.map(geometry, "opacity:0.5;fill:rgb(0,128,0);stroke:rgb(0,64,0);stroke-width:2");
        }
        else
        {
            // TODO: clip input points/linestring
            mapper.map(geometry, "opacity:0.5;stroke:rgb(0,128,0);stroke-width:10");
        }

        {
            // Map buffer in yellow (inflate) and with orange-dots (deflate)
            std::string style = negative
                ? "opacity:0.4;fill:rgb(255,255,192);stroke:rgb(255,128,0);stroke-width:3"
                : "opacity:0.4;fill:rgb(255,255,128);stroke:rgb(0,0,0);stroke-width:3";

            mapper.map(buffered, style);
        }
    }

    template <typename Mapper, typename Geometry, typename GeometryBuffer>
    void map_io_zoomed(Mapper& mapper, Geometry const& geometry,
            GeometryBuffer const& buffered, bool negative, bool areal)
    {
        // Map input geometry in green
        if (areal)
        {
            // Assuming input is areal
            GeometryBuffer clipped;
// TODO: the next line does NOT compile for multi-point, TODO: implement that line
//            bg::intersection(geometry, m_alternate_box, clipped);
            mapper.map(clipped, "opacity:0.5;fill:rgb(0,128,0);stroke:rgb(0,64,0);stroke-width:2");
        }
        else
        {
            // TODO: clip input (multi)point/linestring
            mapper.map(geometry, "opacity:0.5;stroke:rgb(0,128,0);stroke-width:10");
        }

        {
            // Map buffer in yellow (inflate) and with orange-dots (deflate)
            std::string style = negative
                ? "opacity:0.4;fill:rgb(255,255,192);stroke:rgb(255,128,0);stroke-width:3"
                : "opacity:0.4;fill:rgb(255,255,128);stroke:rgb(0,0,0);stroke-width:3";

            try
            {
                // Clip output multi-polygon with box
                GeometryBuffer clipped;
                bg::intersection(buffered, m_alternate_box, clipped);
                mapper.map(clipped, style);
            }
            catch (...)
            {
                std::cout << "Error for buffered output " << m_casename << std::endl;
            }
        }
    }

    bool m_zoom{false};
    bg::model::box<Point> m_alternate_box;
    std::string m_casename;
};


#endif
