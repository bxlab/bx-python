#!/usr/bin/env python2.4

from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib.numerix import arange, sin, pi
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import Affine, Bbox, Value, Point, get_bbox_transform, unit_bbox, identity_transform

# Use GTKAgg directly
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import gobject, gtk
import sys

from bx.tracks import *
from bx.tracks.align import *

import bx.align.maf

def build_tc():

    tc = TrackManager()

    # Build the maf_track
    reader = bx.align.maf.Reader( sys.stdin )
    block = reader.next()
    lo, hi = 0, block.text_size
    t = AlignmentTrack( "Test Alignment", block )
    
    u = AlignmentTrackUnderlay( "Test underlay" )
    u.add_interval( 'rfbat.1', 20, 40, 0.5 )
    u.add_interval( 'rfbat.1', 50, 60, 1.0 )
    u.set_patch_attributes( facecolor='red', linewidth=0 )
    
    t.add_underlay( u )

    tc.set_range( lo, hi )
    tc.add_track( t )
    
    return tc

def set_axes_position( ax, left, bottom, right, top ):
    left, bottom, right, top = map( Value, ( left, bottom, right, top ) )
    ax.left = left
    ax.bottom = bottom
    ax.right = right
    ax.top = top
    ax.bbox = Bbox( Point( left, bottom ), Point( right, top ) )
    ax._set_lim_and_transforms()

def main():
    
    win = gtk.Window()
    win.connect("destroy", lambda x: gtk.main_quit())
    win.set_default_size(600,400)
    win.set_title("Embedding in GTK")

    # Use a Vbox to stack up the widgets in the window
    hpaned = gtk.HPaned()
    win.add( hpaned )

    # Build the actual figure
    tc = build_tc()
    fig = tc.build_figure()

    # Wrap the figure in a GTK widget
    # TODO: Make size of figure depend on size of tracks
    canvas = FigureCanvas( fig )
    tc.canvas = canvas
    
    canvas.set_size_request( fig.ur.x().get(), fig.ur.y().get() )

    # Wrap canvas in an 'Alignment' to prevent it from stretching
    alignment = gtk.Alignment( 0.5, 0.5  )
    alignment.add( canvas )
    
    # Wrap the alignment in a scroll area (so it scrolls when window is smaller than figure)
    scroll = gtk.ScrolledWindow()
    scroll.add_with_viewport( alignment )
    scroll.set_shadow_type( gtk.SHADOW_NONE )
    scroll.get_child().set_shadow_type( gtk.SHADOW_NONE )
    hpaned.add2( scroll )
    
    # Build a tree model
    model = gtk.TreeStore( gobject.TYPE_STRING, gobject.TYPE_PYOBJECT )
    fill_tree_model( model, None, tc.tracks )

    treeview = gtk.TreeView( model )   
    renderer = gtk.CellRendererText()
    column = gtk.TreeViewColumn( "Track", renderer, text=0 )
    treeview.append_column(column)
    
    selection = treeview.get_selection()
    selection.set_mode(gtk.SELECTION_BROWSE)
    
    def treeview_callback( treeview, path, col ):
        model = treeview.get_model()
        iter = model.get_iter( path )
        track = model.get_value( iter, 1 )
        track.do_dialog( tc )
    
    treeview.connect( 'row-activated', treeview_callback)
    treeview.set_size_request(200, -1)

    # Create scrollbars around the view.
    scrolled = gtk.ScrolledWindow()
    scrolled.add(treeview)
    hpaned.add1( scrolled )

    # toolbar = NavigationToolbar(canvas, win)
    # vbox.pack_start(toolbar, False, False)

    win.show_all()
    
    gtk.main()

def fill_tree_model( model, parent, tracks ):
    for t in tracks:
        iter = model.append( parent )
        model.set( iter, 0, t.get_name() )
        model.set( iter, 1, t )
        children = t.get_children()
        if children:
            fill_tree_model( model, iter, children )

if __name__ == "__main__":
    main()
