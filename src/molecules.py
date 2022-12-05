import MDAnalysis as mda
import numpy as np
import sys


class Molecules():

	protein_key = "protein and name CA"
	nucleic_key = "nucleic and name DB"
	

	def __init__(self, topo: str, traj: str):
		self.univ = mda.Universe(topo, traj, trajectory = True)
		self.ca_indices = self.univ.select_atoms(self.protein_key).indices
		self.db_indices = self.univ.select_atoms(self.nucleic_key).indices
		dna_resnames = self.univ.select_atoms(self.nucleic_key).segments.resnames[0] if len(self.db_indices) > 0 else np.array([])
		pro_resnames = self.univ.select_atoms(self.protein_key).resnames
		self.resnames = np.concatenate([dna_resnames, pro_resnames])
		self.framenum = len(self.univ.trajectory)


	def frame(self, iframe: int):
		ts = self.univ.trajectory[iframe]
		ca_pos = ts.positions[self.ca_indices]
		db_pos = ts.positions[self.db_indices]
		reslen = int(len(db_pos) / 2)

		if len(db_pos) % 2 != 0:
			print("Error: The number of dna residues must be even.")
			sys.exit()

		db_coms = 0.5 * (db_pos[:reslen] + db_pos[reslen:][::-1])

		return np.concatenate([db_coms, ca_pos])


