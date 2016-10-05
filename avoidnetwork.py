import requests
import inspect
import logging
logging.basicConfig(level=logging.INFO) #, format='%(levelname)s: %(message)s')

original_get = requests.get

def plotly_no_get(*args, **kwargs):
	one_frame_up = inspect.stack()[1]
	if one_frame_up[3] == 'get_graph_reference':
		raise requests.exceptions.RequestException
	return original_get(*args, **kwargs)

requests.get = plotly_no_get

import plotly
