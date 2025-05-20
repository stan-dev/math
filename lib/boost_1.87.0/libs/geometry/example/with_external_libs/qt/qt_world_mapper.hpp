// Boost.Geometry
//
// Copyright (c) 2007-2024 Barend Gehrels, Amsterdam, the Netherlands.
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Qt World Mapper Example

#ifndef QT_WORLD_MAPPER_H
#define QT_WORLD_MAPPER_H

#include <QWidget>

#include <boost/geometry/geometry.hpp>

#include <boost/geometry/geometries/geometries.hpp>

using point_2d = boost::geometry::model::d2::point_xy<double>;
using country_type = boost::geometry::model::multi_polygon
    <
        boost::geometry::model::polygon<point_2d>
    >;

class qt_world_mapper : public QWidget
{
    Q_OBJECT

public:
    qt_world_mapper(std::vector<country_type> const& countries, boost::geometry::model::box<point_2d> const& box, QWidget *parent = nullptr);

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    using map_transformer_type = boost::geometry::strategy::transform::map_transformer
        <
            double, 2, 2,
            true, true
        >;

     std::vector<country_type> m_countries;
     boost::geometry::model::box<point_2d> m_box;
};


#endif
