from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.models import HoverTool
from bokeh.io import show, output_notebook
from numpy import linspace, array
from scipy.optimize import curve_fit
output_notebook()

def f_scaling(x, y):
    def amdahl(z, a, b): return a + b/z
    par, cov = curve_fit(amdahl, x, y)
    def f(z): return amdahl(z, *par)
    return f, par

bw_times = {
            "awp" : 1051,
            "cactus" : 4800,
            "milc" : 7916,
            "namd" : 242.2,
            "nwchem" : 24160,
            "ppm" : 7790,
            "psdns" : 1538,
            "qmcpack" : 1832,
            "rmg" : 7310,
            "vpic" : 4218,
            "wrf" : 3289
            }
nodes = {
        "awp" : [128, 683, 1067, 1195, 1366],
        "cactus" : [1000, 1200, 1600, 1700],        
        "milc" : [128, 256, 512, 768, 1296],
        "namd" : [512, 1024, 1312, 1600],
        "nwchem" : [800, 1024],
        "ppm" : [704, 1408],
        "psdns" : [460, 820, 1640],
        "qmcpack" : [400, 625, 800, 900, 1000, 1200, 1300, 1600],
        "rmg" : [576, 864, 1152, 1576],
        "vpic" : [512, 572, 768, 1024, 1536],
        "wrf" : [400, 800, 1600, 1680]
        }
times = {
        "awp" : [3677, 699, 502, 402, 335],
        "cactus" : [2258, 2001, 1939, 1879],
        "milc" : [9749, 4722, 2486, 2258, 1364],
        "namd" : [270, 136.5, 112.3, 92.02],
        "nwchem" : [22623, 16884],
        "ppm" : [13410, 7084],
        "psdns" : [2684, 1496, 1009],
        "qmcpack" : [4080, 2726, 2141, 1954, 1768, 1498, 1410, 1334],
        "rmg" : [2640, 2172, 2080, 1880],
        "vpic" : [9183, 8546, 6232, 4573, 3239],
        "wrf" : [6434, 3446, 1768, 1708]
        }
hover = HoverTool(tooltips = [("nodes", "@x"), ("runtime(s)", "@y")])
plots = []
for spp in nodes.keys():
    f, par = f_scaling(nodes[spp], times[spp])
    x = linspace(min(nodes[spp]), 8000)
    p = figure(title = spp, height = 200,width = 300, tools = ["pan,wheel_zoom,box_zoom,reset", hover])
    p.line(x, f(x))
    p.line(x, len(x)*[bw_times[spp]/3.0], color = "red", legend = "T_BW/3")
    p.circle(nodes[spp], times[spp], color = "red")
    bw_times_nodes = par[1]/(bw_times[spp]/3.0-par[0])
    p.circle([bw_times_nodes], [bw_times[spp]/3.0], color = "purple", 
                 legend = "#Nodes " + str(int(bw_times_nodes+0.5)))
    p.legend.glyph_width = 10
    p.legend.glyph_height = 10
    p.legend.label_text_font_size = "8pt"
    plots += [p]
show(gridplot(*plots, ncols=3))
