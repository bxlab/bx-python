from bx.tracks import *

from matplotlib.patches import *
from matplotlib.transforms import *
from matplotlib.colors import colorConverter
from bx.matplotlib.extras import *

from itertools import *

class AlignmentTrack( Track ):
    """
    Draws an alignment block
    """
    def __init__( self, name, block ):
        Track.__init__( self, name )
        self.block = block
        self.underlays = []
        
    def add_underlay( self, underlay ):
        self.underlays.append( underlay )
        
    def get_children( self ):
        return self.underlays
            
    def get_height( self, row_height ):
        return len( self.block.components ) * row_height
        
    def draw( self, ax, xmin, xmax ):
        """
        Draw the MAF into an axes
        """
        # Input bbox is the range of data cordinates
        boxin = Bbox( Point( Value( xmin ), Value( 0 ) ), Point( Value( xmax ), Value( 1 ) ) )
        # We will draw a row for each component
        nrows = len( self.block.components )
        # This is the row height in terms of device cordinates
        row_height = ( ax.bbox.ur().y() - ax.bbox.ll().y() ) / Value( nrows )
        # Build a transformation for each row and then draw it in
        ticklocs, ticklabels = [], []
        for i, comp in enumerate( self.block.components ):
            # Out box is a 1/nrows chunk of the total area
            boxout = Bbox( Point(ax.bbox.ll().x(), ax.bbox.ll().y() + Value( nrows - i - 1 ) * row_height ),
                           Point(ax.bbox.ur().x(), ax.bbox.ll().y() + Value( nrows - i ) * row_height ) )
            # Make a transformation between the two boxes
            trans = get_bbox_transform( boxin, boxout)
            # Todo, draw any 'underlays'
            for underlay in self.underlays:
                underlay.draw( comp.src, ax, trans )
            # Draw each character of the alignment onto the axes
            cc = CharacterCollection( comp.text, 
                                      [ j + 0.5 for j in range( len( comp.text ) ) ], 
                                      0.5, 
                                      horizontalalignment='center', 
                                      verticalalignment='center' ) 
            cc.fontproperties.set_family( 'monospace' )
            cc.set_transform( trans )
            ax.add_artist( cc )
            ## for j, ch in enumerate( comp.text ):
            ##   ax.text( j + 0.5, 0.5, ch,
            ##             verticalalignment='center', 
            ##             horizontalalignment='center',
            ##             transform=trans,
            ##             fontdict=dict(family='monospace') )
            # Add a label
            offset = ( nrows - i ) / (nrows+0.0) - ( 1 / ( 2.0 * nrows ) )
            ticklocs.append( offset )
            ticklabels.append( comp.src )
        # Set labels
        ax.yaxis.set_ticks( ticklocs )
        ax.yaxis.set_ticklabels( ticklabels )
        for t in ax.yaxis.get_major_ticks(): t.tick1On = t.tick2On = False
        ax.xaxis.set_ticks( [] )
            
class AlignmentTrackUnderlay( Track ):
    def __init__( self, name, threshold=0.0 ):
        Track.__init__( self, name )
        self.intervals = {}
        self.threshold = threshold
        self.patches = []
        self.patch_scores = []
    def set_patch_attributes( self, **kwargs ):
        self.patch_attributes = kwargs        
    def add_interval( self, src, start_col, end_col, score=1.0 ):
        try:
            self.intervals[src].append( ( start_col, end_col, score ) )
        except KeyError:
            self.intervals[src] = [ ( start_col, end_col, score ) ]
    def draw( self, src, axes, transform ):
        if src in self.intervals:
            for start_col, end_col, score in self.intervals[src]:
                r = Rectangle( ( start_col, 0 ), end_col - start_col, 1, **self.patch_attributes )
                r.set_transform( transform )
                axes.add_patch( r )
                r.set_visible( score >= self.threshold )
                self.patches.append( r )
                self.patch_scores.append( score )
                
    def set_threshold( self, threshold ):
        self.threshold = threshold
        for patch, score in izip( self.patches, self.patch_scores ):
            patch.set_visible( score >= self.threshold )

    def do_dialog( self, tc ):
        dlg = AlignmentTrackUnderlayDialog( tc, self )
        dlg.show()
        dlg.run()
        
class AlignmentTrackUnderlayDialog(gtk.Dialog):
    def __init__(self, tc, underlay ):
        gtk.Dialog.__init__(self, 'Edit Underlay')

        self.tc = tc
        self.underlay = underlay
        
        self.rgb = colorConverter.to_rgb(self.underlay.patches[0].get_facecolor())
        def set_color(button):
            rgb = get_color(self.rgb)
            if rgb is not None:
                self.rgb = rgb
        button = gtk.Button(stock=gtk.STOCK_SELECT_COLOR)
        button.show()
        button.connect('clicked', set_color)
        
        self.vbox.pack_start( button )

        self.slider_adjustment = gtk.Adjustment( underlay.threshold, 0.0, 1.0 )
        slider = gtk.HScale( self.slider_adjustment )
        slider.show()
        
        self.vbox.pack_start( slider )

        self.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        self.add_button(gtk.STOCK_APPLY, gtk.RESPONSE_APPLY)
        self.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)

    def update(self):
        for patch in self.underlay.patches:
            patch.set_facecolor(self.rgb)
        self.underlay.set_threshold( float( self.slider_adjustment.get_value() ) )
        self.tc.redraw()

    def run(self):
        while 1:
            response = gtk.Dialog.run(self)
            if response==gtk.RESPONSE_APPLY:
                self.update()
            elif response==gtk.RESPONSE_OK:
                self.update()
                break
            elif response==gtk.RESPONSE_CANCEL:
                break
        self.destroy()        
        
# This utility method for a color picker dialog should be shared
         
def get_color(rgb):
    def rgb_to_gdk_color(rgb):
        r,g,b = rgb
        color = gtk.gdk.Color(int(r*65535), int(g*65535), int(b*65535))
        return color

    def gdk_color_to_rgb(color):
        return color.red/65535.0, color.green/65535.0, color.blue/65535.0

    dialog = gtk.ColorSelectionDialog('Choose color')

    colorsel = dialog.colorsel
    color = rgb_to_gdk_color(rgb)
    colorsel.set_previous_color(color)
    colorsel.set_current_color(color)
    colorsel.set_has_palette(True)

    response = dialog.run()

    if response == gtk.RESPONSE_OK:
        rgb = gdk_color_to_rgb(colorsel.get_current_color())
    else:
        rgb = None
    dialog.destroy()
    return rgb  
        
