"""alignment visualization

Took inspiration from Aliview and web 
https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner.
Using Bokeh to vsualiza the alignmnts.

This file can also be imported as a module and contains the 
following main functions:

    * read_fasta - returns list of alignments
    * view_alignment - returns Bokeh image
    * view_alignment_with_time - returns Bokeh image
    
"""


import os, io, random
import string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import operator
import matplotlib.lines as mlines
import pylab as pl


from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, Legend, LegendItem
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.colors import RGB
from bokeh.io import output_file, show, export_png, export_svgs


from sklearn.preprocessing import minmax_scale

import warnings
warnings.filterwarnings('ignore')

def read_fasta(filename):
    """read faseta file

    Parameters
    ----------
    filename : str
     
    Returns
    -------
    alignments: list
        returns list of alignments
  
    """
    read = []
    alignments = {}
    max_len = 0
    f = open(filename, "r")
    for line in f:
        read.append(line)
        length = len(line)
        if length > max_len and '>' not in line:
            max_len = length
    f.close() 
    
    for i in range(0,len(read),2):
        num = max_len - len(read[i+1])
        seq = read[i+1][:-1] + num*'-'
        alignments[read[i][1:-1]] = seq
        
    return alignments

def get_random_color():
    r = random.randint(0,255)
    g = random.randint(0,255)
    b = random.randint(0,255)
    rgb = RGB(r,g,b,0.7)
    return rgb


def get_colors(seqs):
    """color encoder

    Parameters
    ----------
    seqs : str
        sequence of activities
        
    Returns
    -------
    colors: list
        returns list of corresponding colors
  
    """
    text = [i for s in list(seqs) for i in s]
    activities =['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '1', '2', '3','-']
    #color_encoder = {}
    #for key in activities:
        #color_encoder[key] = get_random_color()
    #clrs =  {'A':RGB(255,0,0,0.5),'T':'green','G':'orange','C':'blue','-':'white'}
    color_encoder ={'a': RGB(103, 2, 210,0.9),
 'b': RGB(103, 2, 210, 0.9),
 'c': RGB(103, 2, 210, 0.9),
 'd': RGB(103, 2, 210, 0.9),
 'e': RGB(103, 2, 210, 0.9),
 
 'f': RGB(4, 49, 210, 0.9),
 'g': RGB(4, 49, 210, 0.9),
 'h': RGB(4, 49, 210, 0.9),
 'i': RGB(4, 49, 210, 0.9),
 'j': RGB(4, 49, 210, 0.9),
 
 'k': RGB(1,141,18, 0.9),
 'l': RGB(1,141,18, 0.9),
 'm': RGB(1,141,18, 0.9),
 'n': RGB(1,141,18, 0.9),
 
 'o': RGB(228,112, 5,  0.9),
 'p': RGB(228,112, 5,  0.9),
 'q': RGB(228,112, 5,  0.9),
 
 'r': RGB(243, 31, 7,0.9),
 's': RGB(243, 31, 7,0.9),
 't': RGB(243, 31, 7,0.9),
 'u': RGB(243, 31, 7,0.9),
 'v': RGB(243, 31, 7,0.9),
 'w': RGB(243, 31, 7,0.9),
 'x': RGB(243, 31, 7,0.9),
 
 'y': RGB(243, 219,5, 0.9),
 'z': RGB(243, 219,5, 0.9),
 '0': RGB(243, 219,5, 0.9),
 '1': RGB(243, 219,5, 0.9),
 '2': RGB(243, 219,5, 0.9),
 
 '-': RGB(256, 256, 256, 1)}
    colors = [color_encoder[i] for i in text]
    return colors


def view_alignment(seqs,ids,plot_width=1000,fontsize="14pt",text_font_size="14pt",height_adjust=50):
    """Bokeh sequence alignment view

    Parameters
    ----------
    fasta : str
        filename of a fasta file
        
    Returns
    -------
    p1: Bokeh image
  
    """
    
    #aln = read_fasta(fasta)
    
    #seqs =[aln[key] for key in (aln)]
    #ids = [key for key in (aln)]
    max_len = len(max(seqs, key = len) )
    for i in range(len(seqs)):
        num = max_len - len(seqs[i])
        seqs[i] = seqs[i] + num*'-'
    #print(seqs)    
    ids = [str(i) for i in ids]
    
    row_number = len(seqs)
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)   
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    DataSource = (gx, gy, recty, text, colors)
    p1 = plot_ColumnDataSource(DataSource, row_number,ids, plot_width=plot_width,fontsize=fontsize,text_font_size=text_font_size,height_adjust=height_adjust)
    return p1
    


def view_alignment_with_time(seqs,ids,time_array,time_range=20):
    """Bokeh alignment visualization with time information 

    Parameters
    ----------
    fasta : str
        filename of a fasta file
    time_array : matrix
        time information 
    time_range : int
        length of time bar allocate to each cell
        
    Returns
    -------
    p1: Bokeh image
  
    """
   
    #aln = read_fasta(fasta) 
    time = np.array(time_array)/ time_range
    time = np.round(time)
    column_max = np.max(time, axis=0)
    color_bar = np.tile(column_max,(time.shape[0],1))
    color_bar = color_bar.flatten()
    time = time.flatten()
    ###modify text and color 
    #seqs =[aln[key] for key in (aln)]
    #ids = [key for key in (aln)]
    max_len = len(max(seqs, key = len) )
    for i in range(len(seqs)):
        num = max_len - len(seqs[i])
        seqs[i] = seqs[i] + num*'-'
    #print(seqs)    
    ids = [str(i) for i in ids]
    
    
    row_number = len(seqs)
    text = [i for s in list(seqs) for i in s]
    time_text = time_array.flatten()
    activities =['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '1', '2', '3']
    color_encoder ={'a': RGB(71, 81, 49, 0.7),
 'b': RGB(64, 114, 231, 0.7),
 'c': RGB(143, 5, 33, 0.7),
 'd': RGB(224, 238, 41, 0.7),
 'e': RGB(12, 206, 241, 0.7),
 'f': RGB(186, 98, 123, 0.7),
 'g': RGB(8, 7, 119, 0.7),
 'h': RGB(15, 33, 246, 0.7),
 'i': RGB(218, 97, 85, 0.7),
 'j': RGB(114, 122, 243, 0.7),
 'k': RGB(98, 126, 93, 0.7),
 'l': RGB(102, 10, 199, 0.7),
 'm': RGB(127, 173, 55, 0.7),
 'n': RGB(61, 113, 185, 0.7),
 'o': RGB(231, 79, 207, 0.7),
 'p': RGB(138, 243, 66, 0.7),
 'q': RGB(144, 251, 237, 0.7),
 'r': RGB(116, 248, 47, 0.7),
 's': RGB(170, 226, 35, 0.7),
 't': RGB(217, 179, 0, 0.7),
 'u': RGB(156, 40, 186, 0.7),
 'v': RGB(106, 223, 57, 0.7),
 'w': RGB(229, 171, 84, 0.7),
 'x': RGB(190, 85, 122, 0.7),
 'y': RGB(48, 168, 136, 0.7),
 'z': RGB(17, 215, 195, 0.7),
 '1': RGB(75, 160, 254, 0.7),
 '2': RGB(2, 91, 121, 0.7),
 '3': RGB(168, 105, 224, 0.7),
 '-': RGB(256, 256, 256, 0.7)}
    #colors = [color_encoder[i] for i in text]
    assert len(color_bar) == len(text) , "The number of sequence is not equal to thr number of text !"
    assert len(color_bar) == len(time) , "The number of sequence is not equal to thr number of time !"
    
    
    colors = []
    white = color_encoder['-']
    for i in range(len(text)):
        color = color_encoder[text[i]]
        if color_bar[i] < 5:
            color_bar[i] = 5
            #print(color_bar[i],time[i])
        if (time[i] > 0) and (time[i]  < color_bar[i]):
            
            for n in range(int(time[i])):
                colors.append(color)
            for n in range(int(color_bar[i] - time[i] )):
                colors.append(white)
        elif time[i] == 0:
            colors.append(color)
            for n in range(int(color_bar[i])-1):
                colors.append(white)
        elif time[i]  == color_bar[i]:
            for n in range(int(time[i])):
                colors.append(color)
        else:
            print(" exceptions here !"    )
    text_modified = []

    for a in range(len(text)):
        text_modified.append(text[a])
        #print(color_bar[a])
        if text[a] == '-':
            length = 1
        else:
            
            text_modified.append(' ')
            text_modified.append(str(time_text[a]))
            text_modified.append(' ')
            text_modified.append('s')
            length = 5
        for i in range(int(color_bar[a] - length)):
            text_modified.append(' ')
    #print(len(text_modified))  
    
    
    N = len(text_modified)/len(seqs)  
    S = len(seqs)    
    width = .4
    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    
    
    #now we can create the ColumnDataSource with all the arrays
    
    DataSource = (gx, gy, recty, text_modified, colors)
    p1 = plot_ColumnDataSource(DataSource, row_number, ids, plot_width=5000)
    return p1
    
def plot_ColumnDataSource(DataSource, row_number, ids, plot_width=1000, fontsize="14pt",text_font_size="14pt",height_adjust=50):
    """Bokeh alignment visualization with time information 

    Parameters
    ----------
    DataSource : tuple
        computed columnDataSource including gx, gy, recty, text_modified and colors
    row_number : int
        number of sequences 
    ids : str
        index, name or information of a sequence
    plot_width : int 
        width of figure
    fontsize : str
        fontsize of figure
    Returns
    -------
    p: Bokeh image
  
    """
    
    gx, gy, recty, text_modified, colors = DataSource
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text_modified, colors=colors))
    #print("Check if these numbers are equal? ",len(gx),len(gy),len(text_modified),len(colors))
    assert len(gx) == len(gy) and len(gx) == len(text_modified) and len(gx) == len(colors), 'error, These Column Data Sources are note equal'
    
    N = len(text_modified)/row_number 
    S = row_number 
    
    
    plot_height = row_number*15+height_adjust
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"
    
    #entire sequence view (no text, with zoom)
    p = figure(title=None, plot_width= plot_width, plot_height=50,
            x_range=x_range, y_range=(0,S), tools=tools,
            min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  
    
    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
            x_range=x_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)    
    
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=text_font_size,text_font_style = "bold")
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)
    
    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.major_label_text_font_style = 'bold'  
    p1.yaxis.major_label_text_font_size = fontsize   
    p1.xaxis.major_label_text_font_size = fontsize 
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
	
    
    #legend = Legend(items=[("test"   , [rects]),], location=(0, -30))

    #p1.add_layout(legend, 'left')
    #li1 = LegendItem(label='Stage 1', renderers=[p1.renderers[1]])
    #li2 = LegendItem(label='Stage 2', renderers=[p1.renderers[1]])
    #li3 = LegendItem(label='Stage 3', renderers=[p1.renderers[1]])
    #li4 = LegendItem(label='Stage 4', renderers=[p1.renderers[1]])
    #li5 = LegendItem(label='Stage 5', renderers=[p1.renderers[1]])
    #li6 = LegendItem(label='Stage 6', renderers=[p1.renderers[1]])
    #legend1 = Legend(items=[li1, li2, li3,li4, li5, li6 ], location='bottom_right')
    #p1.add_layout(legend1)
    #p1.legend.label_text_font_size = '10pt'
   #print(p1.renderers)
    p = gridplot([[p],[p1]], toolbar_location='below')
    return p1