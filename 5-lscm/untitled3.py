import pywavefront
from pywavefront import visualization

obj = pywavefront.Wavefront('input-face.obj')

visualization.draw(obj)