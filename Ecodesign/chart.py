#code source:https://plotly.com/python/images/
import plotly.graph_objects as go
from PIL import Image

def screen():
	fig = go.Figure()
# Add image
	img_width = 1600
	img_height = 900
	scale_factor = 0.5

	filename = "Figures\Toughness_Price_Chart.png"
    
	fig.add_layout_image(
        	x=0,
        	sizex=img_width,
        	y=0,
        	sizey=img_height,
        	xref="x",
        	yref="y",
        	opacity=1.0,
        	layer="below",
        	source= Image.open(filename)
	)
	fig.update_xaxes(showgrid=False, range=(0, img_width))
	fig.update_yaxes(showgrid=False, scaleanchor='x', range=(img_height, 0))

# Set dragmode and newshape properties; add modebar buttons
	fig.update_layout(
    		dragmode='drawrect',
    		newshape=dict(line_color='black'),
   		title_text='Drag to add annotations - use modebar to change drawing tool')
	fig.show(config={'modeBarButtonsToAdd':['drawline',
                                        'drawopenpath',
                                        'drawclosedpath',
                                        'drawcircle',
                                        'drawrect',
                                        'eraseshape'
                                       ]})

def rank():
	fig = go.Figure()
# Add image
	img_width = 1600
	img_height = 900
	scale_factor = 0.5

	filename = "Figures\Modulus_Density_Chart.png"
    
	fig.add_layout_image(
        	x=0,
        	sizex=img_width,
        	y=0,
        	sizey=img_height,
        	xref="x",
        	yref="y",
        	opacity=1.0,
        	layer="below",
        	source= Image.open(filename)
	)
	fig.update_xaxes(showgrid=False, range=(0, img_width))
	fig.update_yaxes(showgrid=False, scaleanchor='x', range=(img_height, 0))

# Set dragmode and newshape properties; add modebar buttons
	fig.update_layout(
    		dragmode='drawrect',
    		newshape=dict(line_color='black'),
   		title_text='Drag to add annotations - use modebar to change drawing tool')
	fig.show(config={'modeBarButtonsToAdd':['drawline',
                                        'drawopenpath',
                                        'drawclosedpath',
                                        'drawcircle',
                                        'drawrect',
                                        'eraseshape'
                                       ]})