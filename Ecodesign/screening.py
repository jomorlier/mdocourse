#https://ipywidgets.readthedocs.io/en/7.6.2/examples/Widget%20Events.html
import ipywidgets as widgets
import IPython
from IPython.display import display
from IPython.display import Video

def results():
	button = widgets.Button(description="View Results")
	output= widgets.Output()
	
	display(button,output)

	def on_button_clicked(b):
		with output:
			display(Video(filename='screeningRV.mp4',height=600))
	button.on_click(on_button_clicked)