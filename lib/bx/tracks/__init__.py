import gtk

from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib.numerix import arange, sin, pi
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import Affine, Bbox, Value, Point, get_bbox_transform, unit_bbox, identity_transform

# Some definitions related to display sizes
DEFAULT_FONT_SIZE = 14.0 # points, should be determined on the fly from the font
DEFAULT_CHAR_WIDTH = 8.0 # this is just a guess, need to determine it from font metrics
POINTS_PER_INCH = 72.0

# Height of a single 'track', based on character height
TRACK_HEIGHT = DEFAULT_FONT_SIZE / POINTS_PER_INCH
# Spacing between tracks
TRACK_SPACE = TRACK_HEIGHT / 4.0

class Track( object ):
    """
    Abstract base class for tracks
    """
    def __init__( self, name ):
        self.name = name
    def get_name( self ):
        return self.name
    def get_children( self ):
        return []
    def do_dialog( self ):
        pass
    
class TrackManager( object ):
    """
    Manages a collection of tracks. Responsible for building a figure
    and drawing / redrawing each track onto it.
    """
    
    def __init__( self ):
        self.tracks = []
        
    def set_range( self, start_pos, stop_pos ):
        self.start_pos = start_pos
        self.stop_pos = stop_pos
        
    def add_track( self, track ):
        self.tracks.append( track )
        
    def build_figure( self ):
    
        total_track_height = sum( [ t.get_height( TRACK_HEIGHT ) for t in self.tracks ] )
        total_track_height += TRACK_SPACE * ( len( self.tracks ) - 1 )
    
        # Build the figure
        margin_l, margin_r, margin_t, margin_b = 2.5, 0.5, 0.5, 0.5
    
        # Width and height determined to fit the alignment
        fig_height = margin_t + margin_b + total_track_height
        fig_width = margin_r + margin_l + DEFAULT_CHAR_WIDTH * self.stop_pos / POINTS_PER_INCH

        # Track height in figure relative coords
        rel_track_height = TRACK_HEIGHT / fig_height
        
        fig = Figure( figsize=(fig_width,fig_height), frameon=True )
    
        top = 1 - ( margin_t / fig_height )
    
        for track in self.tracks:
    
            track_height = track.get_height( rel_track_height )
    
            ax_left = margin_l / fig_width
            ax_width = ( fig_width - margin_l - margin_r ) / fig_width
            ax_height = track_height
            ax_bottom = top - ax_height
            
            ax = fig.add_axes( ( ax_left, ax_bottom, ax_width, ax_height ), frameon=False )
            track.draw( ax, self.start_pos, self.stop_pos )
            ax.set_xlim( ( self.start_pos, self.stop_pos ) )
    
            top -= ( ax_height + TRACK_SPACE )
    
        self.fig = fig
        return fig
        
    def redraw( self ):
        self.canvas.draw()
