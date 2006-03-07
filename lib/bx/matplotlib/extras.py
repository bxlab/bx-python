from matplotlib.artist import Artist
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties

class CharacterCollection(Artist):
    """
    A whole bunch of characters with the same font, alignment, color,
    et cetera, but with different positions
    """

    def __init__( self, c, x, y, 
                  verticalalignment='bottom',
                  horizontalalignment='left',
                  color=None, 
                  fontproperties=None ):
        """
        c: list of characters to draw
        x: list of x positions OR single shared x position
        y: list of y positions OR single shared y position
        """
        Artist.__init__( self )
        self.c = c
        self.x = x
        self.y = y
        self.color = color or rcParams['text.color']
        self.fontproperties = fontproperties or FontProperties()
        self.angle = 0
        self.verticalalignment = verticalalignment
        self.horizontalalignment = horizontalalignment
            
        
    def draw(self, renderer):
        if not self.get_visible(): return
        renderer.open_group('charcollection')
        
        gc = renderer.new_gc()
        gc.set_foreground(self.color)
        gc.set_alpha(self._alpha)

        is_x_seq = isinstance( self.x, list )
        is_y_seq = isinstance( self.y, list )

        # This gives us the T-box for vertical alignment 
        tw, th = renderer.get_text_width_height( 'T', self.fontproperties, False )

        # Cache the metrics for a given string, speeds up things like
        # DNA sequence rendering (~9 characters, over and over)
        cache = {}
        # Draw each character in the sequence
        for i in range( len( self.c ) ):
            # Get c, x, and y
            c = self.c[i]
            if is_x_seq: x = self.x[i]
            else: x = self.x
            if is_y_seq: y = self.y[i]
            else: y = self.y
            # Get ascent and descent
            if c in cache: 
                w, h, a, d = cache[c]
            else:
                w, h = renderer.get_text_width_height( c, self.fontproperties, False )
                d = renderer.get_text_descent( c, self.fontproperties, False )
                a = th - ( h - d )
                cache[c] = w, h, a, d
            # Determine transformed positions
            tx, ty = self._transform.xy_tup( ( x, y ) )
            ## print "Transformed", c, tx, ty, d
            ## renderer.draw_rectangle( gc, None, tx, ty, w, h )
            # Offsets for horizontal alignment
            if self.horizontalalignment=='center':  offsetx = tx - (w/2.0)
            elif self.horizontalalignment=='right': offsetx = tx - (w)
            else: offsetx = tx
            # Offsets for vertical alignment
            if self.verticalalignment=='center': offsety = ty - d - th/2.0 # ((th-d+a)/2.0) - d
            elif self.verticalalignment=='top': offsety  = ty - ( h + a )
            else: offsety = ty - d 
            # I really do not understand why the renderer doesn't just handle this itself. 
            if renderer.flipy():
                canvasw, canvash = renderer.get_canvas_width_height()
                offsety = canvash-offsety
            # Draw it
            renderer.draw_text(gc, offsetx, offsety, c, self.fontproperties, self.angle, False )
        renderer.close_group('charcollection')