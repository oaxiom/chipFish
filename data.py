"""
data.py


"""

valid_track_draw_types = frozenset(["graph", "graph_split_strand", "bar", "spot", "kde_graph", "genome"])

valid_move_modes = frozenset(["left", "right", "zoomin", "zoomout"])

colour_lookup_name = {"red": (1, 0, 0),
    "green": (0, 1, 0),
    "blue": (0, 0, 1),
    "orange": (1,0.5,0)
    }