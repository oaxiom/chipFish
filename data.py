"""
data.py


"""

from matplotlib import colors

valid_track_draw_types = frozenset(["graph", "graph_split_strand", "bar", "spot", "kde_graph", "genome", "repeats", 'splice', 'genome_sql'])

valid_move_modes = frozenset(["left", "right", "zoomin", "zoomout"])

def hex_to_rgb(h):
    return tuple(int(h.lstrip('#')[i:i+2], 16)/256 for i in (0, 2 ,4))

colour_lookup_name = colors.ColorConverter.colors
for k in  colors.ColorConverter.colors:
    if '#' in  colors.ColorConverter.colors[k]:
        colour_lookup_name[k] = hex_ = hex_to_rgb(colors.ColorConverter.colors[k])

