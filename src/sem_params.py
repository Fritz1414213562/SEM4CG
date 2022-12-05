from typing import Tuple

class SEMParams():

	grid_params: Tuple[float] = ...
	box_params:  Tuple[float] = ...
	pore_params: Tuple[float] = ...
	cutoff: float = ...
	volt_bottom: float = ...
	volt_top: float = ...
	conductivity_path = "../dat/conductivity_table_amber14sb.h5"

	def __init__():
		pass