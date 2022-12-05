import numpy as np
import sys
from typing import Tuple
#from conductivity import Residue2Conductivity, ConductivityOf
from conductivity import Residue2Conductivity
from sem_params import SEMParams


class GridSpace():

	MAXDIST = 100000

	def __init__(self, box_size: Tuple[float], grid_width: Tuple[float]):
		self.size = box_size
		self.width = grid_width
		self.grids = np.array([])
		self.xs = np.arange(-0.5 * self.size[0], 0.5 * self.size[0] + self.width[0], self.width[0])
		self.ys = np.arange(-0.5 * self.size[1], 0.5 * self.size[1] + self.width[1], self.width[1])
		self.zs = np.arange(-0.5 * self.size[2], 0.5 * self.size[2] + self.width[2], self.width[2])


	def set_pore_params(self, pore_length, upper_width, lower_width):

		X, Y, Z = np.meshgrid(self.xs, self.ys, self.zs, indexing = 'ij')

		minor_rad = 0.25 / (upper_width - lower_width) * (pore_length ** 2 + (upper_width - lower_width) ** 2)
		major_rad = minor_rad + 0.5 * lower_width

		if SEMParams.log_flag:
			print("	Nanopore parameters")
			print("	Pore length: {}".format(SEMParams.pore_params[0]))
			print("	Major Radius: {}".format(major_rad))
			print("	Minor Radius: {}".format(minor_rad))

		R = np.sqrt(X ** 2 + Y ** 2)
		self.grids = np.zeros_like(X)
		self.grids[np.where((R - major_rad) ** 2 + Z ** 2 > minor_rad ** 2)] = 1
		self.grids[np.where(R > major_rad)] = 0
		self.grids[np.where(Z < - 0.5 * pore_length)] = 1
		self.grids[np.where(Z >   0.5 * pore_length)] = 1


	def compute_conductivity(self, cutoff, positions, resnames, bulk):

		if len(positions) > 0:
			dists = np.zeros_like(self.grids, dtype = [("value", float), ("function", object)])
			dists["value"] = self.MAXDIST * self.grids
			dists["function"] = Residue2Conductivity["BULK"]

			posx = positions[:, 0]
			posy = positions[:, 1]
			posz = positions[:, 2]

			xyz_init = -0.5 * self.size
			xyz_fin  =  0.5 * self.size

			xi = np.maximum(np.around((posx - cutoff) / self.width[0]) * self.width[0], xyz_init[0] * np.ones_like(posx))
			xi = np.where(np.abs(xi - posx) < cutoff, np.maximum(xi - self.width[0], xyz_init[0] * np.ones_like(posx)), xi)
			xi_index = np.around((xi - xyz_init[0]) / self.width[0]).astype(int)
			xf = np.minimum(np.around((posx + cutoff) / self.width[0]) * self.width[0], xyz_fin[0] * np.ones_like(posx))
			xf = np.where(np.abs(xf - posx) < cutoff, np.minimum(xf + self.width[0], xyz_fin[0] * np.ones_like(posx)), xf)
			xf_index = np.around((xf - xyz_init[0]) / self.width[0]).astype(int)
			del(xi)
			del(xf)

			yi = np.maximum(np.around((posy - cutoff) / self.width[1]) * self.width[1], xyz_init[1] * np.ones_like(posy))
			yi = np.where(np.abs(yi - posy) < cutoff, np.maximum(yi - self.width[1], xyz_init[1] * np.ones_like(posy)), yi)
			yi_index = np.around((yi - xyz_init[1]) / self.width[1]).astype(int)
			yf = np.minimum(np.around((posy + cutoff) / self.width[1]) * self.width[1], xyz_fin[1] * np.ones_like(posy))
			yf = np.where(np.abs(yf - posy) < cutoff, np.minimum(yf + self.width[1], xyz_fin[1] * np.ones_like(posy)), yf)
			yf_index = np.around((yf - xyz_init[1]) / self.width[1]).astype(int)
			del(yi)
			del(yf)

			zi = np.maximum(np.around((posz - cutoff) / self.width[2]) * self.width[2], xyz_init[2] * np.ones_like(posz))
			zi = np.where(np.abs(zi - posz) < cutoff, np.maximum(zi - self.width[2], xyz_init[2] * np.ones_like(posz)), zi)
			zi_index = np.around((zi - xyz_init[2]) / self.width[2]).astype(int)
			zf = np.minimum(np.around((posz + cutoff) / self.width[2]) * self.width[2], xyz_fin[2] * np.ones_like(posz))
			zf = np.where(np.abs(zf - posz) < cutoff, np.minimum(zf + self.width[2], xyz_fin[2] * np.ones_like(posz)), zf)
			zf_index = np.around((zf - xyz_init[2]) / self.width[2]).astype(int)
			del(zi)
			del(zf)

			retval = bulk * self.grids

			for iatom in range(0, len(positions)):

				surrounds = np.array([\
					[\
						[\
							[x, y, z] for z in self.zs[zi_index[iatom]:zf_index[iatom]]\
						] for y in self.ys[yi_index[iatom]:yf_index[iatom]]\
					] for x in self.xs[xi_index[iatom]:xf_index[iatom]]\
				])

				surrounds = np.sqrt(np.sum(np.square(surrounds - positions[iatom]), axis = 3))
				#surr_indices = np.zeros_like(dists, dtype = bool)
				#surr_indices[xi_index[iatom]:xf_index[iatom], yi_index[iatom]:yf_index[iatom], zi_index[iatom]:zf_index[iatom]] = True
				#dist_section = dists[surr_indices]
				cond_func = Residue2Conductivity[resnames[iatom]][0]
				dist_section = dists[xi_index[iatom]:xf_index[iatom], yi_index[iatom]:yf_index[iatom], zi_index[iatom]:zf_index[iatom]]
				isClose = (surrounds < cutoff) & (surrounds < dist_section["value"])
				dists[xi_index[iatom]:xf_index[iatom], yi_index[iatom]:yf_index[iatom], zi_index[iatom]:zf_index[iatom]]["value"] =\
					np.where(isClose, surrounds, dist_section["value"])
				dists[xi_index[iatom]:xf_index[iatom], yi_index[iatom]:yf_index[iatom], zi_index[iatom]:zf_index[iatom]]["function"] =\
					np.where(isClose, cond_func, dist_section["function"])


			#	for z in dists:
			#		for y in z:
			#			for x in y:
			#				print("{:}".format(x))
			
			retval[(dists["value"] < self.MAXDIST - 1) & (dists["value"] > 0)] = np.vectorize(lambda x: x[1](x[0]))(dists[(dists["value"] < self.MAXDIST - 1) & (dists["value"] > 0)])

			return retval

		return bulk * self.grids
	

