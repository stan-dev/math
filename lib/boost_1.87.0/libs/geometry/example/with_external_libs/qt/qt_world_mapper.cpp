// Boost.Geometry
//
// Copyright (c) 2007-2024 Barend Gehrels, Amsterdam, the Netherlands.
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Qt World Mapper Example

#include "qt_world_mapper.hpp"
#include "common/read_countries.hpp"

#include <QApplication>
#include <QPainter>
#include <QTime>
#include <QTimer>

#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>

// Adapt a QPointF such that it can be handled by Boost.Geometry
BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(QPointF, double, cs::cartesian, x, y, setX, setY)

// Adapt a QPolygonF as well.
// A QPolygonF has no holes (interiors) so it is similar to a Boost.Geometry ring
BOOST_GEOMETRY_REGISTER_RING(QPolygonF)


qt_world_mapper::qt_world_mapper(std::vector<country_type> const& countries, boost::geometry::model::box<point_2d> const& box, QWidget *parent)
    : QWidget(parent)
    , m_countries(countries)
    , m_box(box)
{
    setPalette(QPalette(Qt::blue));
    setAutoFillBackground(true);
}

void qt_world_mapper::paintEvent(QPaintEvent *)
{
    map_transformer_type transformer(m_box, this->width(), this->height());

    QPainter painter(this);
    painter.setBrush(Qt::green);
    painter.setRenderHint(QPainter::Antialiasing);

    for(auto const& country : m_countries)
    {
        for(auto const& polygon : country)
        {
            // This is the essention:
            // Directly transform from a multi_polygon (ring-type) to a QPolygonF
            QPolygonF qring;
            boost::geometry::transform(boost::geometry::exterior_ring(polygon), qring, transformer);
            painter.drawPolygon(qring);
        }
    }
}

int main(int argc, char *argv[])
{
    const std::string filename = argc > 1 ? argv[1] : "../../../data/world.wkt";
    const auto countries = read_countries<country_type>(filename);
    if (countries.empty())
    {
        std::cout << "No countries read" << std::endl;
        return 1;
    }

    const auto box = calculate_envelope<boost::geometry::model::box<point_2d>>(countries);

    QApplication app(argc, argv);
    qt_world_mapper mapper(countries, box);
    mapper.setWindowTitle("Boost.Geometry for Qt - Hello World!");
    mapper.setGeometry(100, 100, 1024, 768);
    mapper.show();
    return app.exec();
}
