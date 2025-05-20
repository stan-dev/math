// Boost.Geometry
//
// Copyright (c) 2010-2024 Barend Gehrels, Amsterdam, the Netherlands.
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// wxWidgets World Mapper example
//
// It will show a basic wxWidgets window, displaying world countries,
// highlighting the country under the mouse, and indicating position
// of the mouse in latitude/longitude and in pixels.

// #define EXAMPLE_WX_USE_GRAPHICS_CONTEXT 1

#include <sstream>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>

#include <wx/wx.h>
#include <wx/math.h>
#include <wx/stockitem.h>
#include <wx/cmdline.h>

#ifdef EXAMPLE_WX_USE_GRAPHICS_CONTEXT
#include "wx/graphics.h"
#include "wx/dcgraph.h"
#endif

#include "common/read_countries.hpp"

using point_2d = boost::geometry::model::d2::point_xy<double>;
using country_type = boost::geometry::model::multi_polygon
    <
        boost::geometry::model::polygon<point_2d>
    >;

// Adapt wxWidgets points to Boost.Geometry points such that they can be used
// in e.g. transformations (see below)
BOOST_GEOMETRY_REGISTER_POINT_2D(wxPoint, int, cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_POINT_2D(wxRealPoint, double, cs::cartesian, x, y)


// ----------------------------------------------------------------------------
class HelloWorldFrame: public wxFrame
{
public:
    HelloWorldFrame(wxFrame *frame, wxString const& title, wxPoint const& pos, wxSize const& size);

    void OnCloseWindow(wxCloseEvent& );
    void OnExit(wxCommandEvent& );

    DECLARE_EVENT_TABLE()
};


// ----------------------------------------------------------------------------
class HelloWorldCanvas: public wxWindow
{
public:
    HelloWorldCanvas(wxFrame *frame, const std::string& filename);

private:
    void DrawCountries(wxDC& dc);
    void DrawCountry(wxDC& dc, country_type const& country);

    void OnPaint(wxPaintEvent& );
    void OnMouseMove(wxMouseEvent&);

    using map_transformer_type = boost::geometry::strategy::transform::map_transformer
        <
            double, 2, 2,
            true, true
        >;

    using inverse_transformer_type = boost::geometry::strategy::transform::inverse_transformer
        <
            double, 2, 2
        >;

    std::shared_ptr<map_transformer_type> m_map_transformer;
    std::shared_ptr<inverse_transformer_type> m_inverse_transformer;

    std::string m_filename;
    boost::geometry::model::box<point_2d> m_box;
    std::vector<country_type> m_countries;
    int m_focus = -1;

    wxBrush m_orange = wxBrush(wxColour(255, 128, 0), wxBRUSHSTYLE_SOLID);
    wxBrush m_blue = wxBrush(wxColour(0, 128, 255), wxBRUSHSTYLE_SOLID);
    wxBrush m_green = wxBrush(wxColour(0, 255, 0), wxBRUSHSTYLE_SOLID);

    wxFrame* m_owner = nullptr;

    DECLARE_EVENT_TABLE()
};



// ----------------------------------------------------------------------------
class HelloWorldApp: public wxApp
{
public:
    bool OnInit()
    {
        const auto filename = wxApp::argc >= 2 ? wxApp::argv[1].ToStdString() : "example/data/world.wkt";

        // Create the main frame window
        HelloWorldFrame *frame = new HelloWorldFrame(NULL, _T("Boost.Geometry for wxWidgets - Hello World!"), wxDefaultPosition, wxSize(1024, 768));

        wxMenu *file_menu = new wxMenu;
        file_menu->Append(wxID_EXIT, wxGetStockLabel(wxID_EXIT));
        wxMenuBar* menuBar = new wxMenuBar;
        menuBar->Append(file_menu, _T("&File"));
        frame->SetMenuBar(menuBar);

        (void) new HelloWorldCanvas(frame, filename);

        // Show the frame
        frame->Show(true);

        return true;
    }
};



// ----------------------------------------------------------------------------
HelloWorldFrame::HelloWorldFrame(wxFrame *frame, wxString const& title, wxPoint const& pos, wxSize const& size)
    : wxFrame(frame, wxID_ANY, title, pos, size, wxDEFAULT_FRAME_STYLE | wxFULL_REPAINT_ON_RESIZE )
{
    CreateStatusBar(2);
}


void HelloWorldFrame::OnExit(wxCommandEvent& )
{
    this->Destroy();
}

void HelloWorldFrame::OnCloseWindow(wxCloseEvent& )
{
    static bool destroyed = false;
    if (! destroyed)
    {
        this->Destroy();
        destroyed = true;
    }
}


// ----------------------------------------------------------------------------
HelloWorldCanvas::HelloWorldCanvas(wxFrame *frame, const std::string& filename)
    : wxWindow(frame, wxID_ANY)
    , m_owner(frame)
    , m_filename(filename)
{
    m_countries = read_countries<country_type>(m_filename);
    m_box = calculate_envelope<boost::geometry::model::box<point_2d>>(m_countries);
}


void HelloWorldCanvas::OnMouseMove(wxMouseEvent &event)
{
    if (m_inverse_transformer)
    {
        // Boiler-plate wxWidgets code
        wxClientDC dc(this);
        PrepareDC(dc);
        m_owner->PrepareDC(dc);

        // Transform the point opn the screen back to Lon/Lat
        point_2d point;
        boost::geometry::transform(event.GetPosition(), point,
                                   *m_inverse_transformer);

        // Determine selected object
        int i = 0;
        int previous_focus = m_focus;
        m_focus = -1;
        for (country_type const& country : m_countries)
        {
            if (boost::geometry::within(point, country))
            {
                m_focus = i;
                break;
            }
            i++;
        }

        // On change:
        if (m_focus != previous_focus)
        {
            // Undraw old focus
            if (previous_focus >= 0)
            {
                dc.SetBrush(m_green);
                DrawCountry(dc, m_countries[previous_focus]);
            }
            // Draw new focus
            if (m_focus >= 0)
            {
                dc.SetBrush(m_orange);
                DrawCountry(dc, m_countries[m_focus]);
            }
        }

        // Create a string and set it in the status text
        std::ostringstream out;
        out << "Position: " << point.x() << ", " << point.y();
        m_owner->SetStatusText(wxString(out.str().c_str(), wxConvUTF8));
    }
}



void HelloWorldCanvas::OnPaint(wxPaintEvent& )
{
    if (m_countries.empty())
    {
        return;
    }

#if defined(EXAMPLE_WX_USE_GRAPHICS_CONTEXT)
    wxPaintDC pdc(this);
    wxGCDC gdc(pdc);
    wxDC& dc = (wxDC&) gdc;
#else
    wxPaintDC dc(this);
#endif

    PrepareDC(dc);

    static bool running = false;
    if (! running)
    {
        running = true;

        // Update the transformers
        wxSize sz = dc.GetSize();
        m_map_transformer.reset(new map_transformer_type(m_box, sz.x, sz.y));
        m_inverse_transformer.reset(new inverse_transformer_type(*m_map_transformer));

        DrawCountries(dc);

        running = false;
    }
}


void HelloWorldCanvas::DrawCountries(wxDC& dc)
{
    namespace bg = boost::geometry;

    dc.SetBackground(m_blue);
    dc.Clear();

    for (country_type const& country :  m_countries)
    {
        dc.SetBrush(m_green);
        DrawCountry(dc, country);
    }
    if (m_focus != -1)
    {
        dc.SetBrush(m_orange);
        DrawCountry(dc, m_countries[m_focus]);
    }
}


void HelloWorldCanvas::DrawCountry(wxDC& dc, country_type const& country)
{
    for (auto const& poly : country)
    {
        // Use only exterior ring, holes are (for the moment) ignored.
        // This would need a holey-polygon compatible wx object

        // Define a Boost.Geometry ring of wxPoints
        // Behind the scenes that is a vector, and a vector has .data(),
        // can be used for the *wxPoint pointer needed for wxWidget DrawPolygon
        boost::geometry::model::ring<wxPoint> ring;
        boost::geometry::transform(boost::geometry::exterior_ring(poly), ring,
                                   *m_map_transformer);
        dc.DrawPolygon(ring.size(), ring.data());
    }
}

// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(HelloWorldFrame, wxFrame)
    EVT_CLOSE(HelloWorldFrame::OnCloseWindow)
    EVT_MENU(wxID_EXIT, HelloWorldFrame::OnExit)
END_EVENT_TABLE()


BEGIN_EVENT_TABLE(HelloWorldCanvas, wxWindow)
    EVT_PAINT(HelloWorldCanvas::OnPaint)
    EVT_MOTION(HelloWorldCanvas::OnMouseMove)
END_EVENT_TABLE()


IMPLEMENT_APP(HelloWorldApp)
