from enum import Enum
from scipy.interpolate import interp1d
import h5py
from const import SEMConst
from sem_params import SEMParams


def bulk_conductivity(_):
	return SEMConst.BULK

#class ConductivityOf(Enum):
#
#	__fin = h5py.File(SEMParams.conductivity_path, 'r')
#	DA   = interp1d(__fin["DA"  + "/r"], __fin["DA"  + "/conductivity"], kind = "cubic"),
#	DT   = interp1d(__fin["DT"  + "/r"], __fin["DT"  + "/conductivity"], kind = "cubic"),
#	DG   = interp1d(__fin["DG"  + "/r"], __fin["DG"  + "/conductivity"], kind = "cubic"),
#	DC   = interp1d(__fin["DC"  + "/r"], __fin["DC"  + "/conductivity"], kind = "cubic"),
#	ALA  = interp1d(__fin["ALA" + "/r"], __fin["ALA" + "/conductivity"], kind = "cubic"),
#	CYS  = interp1d(__fin["CYS" + "/r"], __fin["CYS" + "/conductivity"], kind = "cubic"),
#	ASP  = interp1d(__fin["ASP" + "/r"], __fin["ASP" + "/conductivity"], kind = "cubic"),
#	GLU  = interp1d(__fin["GLU" + "/r"], __fin["GLU" + "/conductivity"], kind = "cubic"),
#	PHE  = interp1d(__fin["PHE" + "/r"], __fin["PHE" + "/conductivity"], kind = "cubic"),
#	GLY  = interp1d(__fin["GLY" + "/r"], __fin["GLY" + "/conductivity"], kind = "cubic"),
#	HIS  = interp1d(__fin["HIS" + "/r"], __fin["HIS" + "/conductivity"], kind = "cubic"),
#	ILE  = interp1d(__fin["ILE" + "/r"], __fin["ILE" + "/conductivity"], kind = "cubic"),
#	LYS  = interp1d(__fin["LYS" + "/r"], __fin["LYS" + "/conductivity"], kind = "cubic"),
#	LEU  = interp1d(__fin["LEU" + "/r"], __fin["LEU" + "/conductivity"], kind = "cubic"),
#	MET  = interp1d(__fin["MET" + "/r"], __fin["MET" + "/conductivity"], kind = "cubic"),
#	ASN  = interp1d(__fin["ASN" + "/r"], __fin["ASN" + "/conductivity"], kind = "cubic"),
#	PRO  = interp1d(__fin["PRO" + "/r"], __fin["PRO" + "/conductivity"], kind = "cubic"),
#	GLN  = interp1d(__fin["GLN" + "/r"], __fin["GLN" + "/conductivity"], kind = "cubic"),
#	ARG  = interp1d(__fin["ARG" + "/r"], __fin["ARG" + "/conductivity"], kind = "cubic"),
#	SER  = interp1d(__fin["SER" + "/r"], __fin["SER" + "/conductivity"], kind = "cubic"),
#	THR  = interp1d(__fin["THR" + "/r"], __fin["THR" + "/conductivity"], kind = "cubic"),
#	VAL  = interp1d(__fin["VAL" + "/r"], __fin["VAL" + "/conductivity"], kind = "cubic"),
#	TRP  = interp1d(__fin["TRP" + "/r"], __fin["TRP" + "/conductivity"], kind = "cubic"),
#	TYR  = interp1d(__fin["TYR" + "/r"], __fin["TYR" + "/conductivity"], kind = "cubic"),
#	BULK = bulk_conductivity
#	__fin.close()

#	@classmethod
#	def set_parameter(cls, parameter_file):
#
#		with h5py.File(parameter_file, 'r') as fin:
#			cls.DA   = interp1d(fin["DA"  + "/r"], fin["DA"  + "/conductivity"], kind = "cubic")
#			cls.DT   = interp1d(fin["DT"  + "/r"], fin["DT"  + "/conductivity"], kind = "cubic")
#			cls.DG   = interp1d(fin["DG"  + "/r"], fin["DG"  + "/conductivity"], kind = "cubic")
#			cls.DC   = interp1d(fin["DC"  + "/r"], fin["DC"  + "/conductivity"], kind = "cubic")
#			cls.ALA  = interp1d(fin["ALA" + "/r"], fin["ALA" + "/conductivity"], kind = "cubic")
#			cls.CYS  = interp1d(fin["CYS" + "/r"], fin["CYS" + "/conductivity"], kind = "cubic")
#			cls.ASP  = interp1d(fin["ASP" + "/r"], fin["ASP" + "/conductivity"], kind = "cubic")
#			cls.GLU  = interp1d(fin["GLU" + "/r"], fin["GLU" + "/conductivity"], kind = "cubic")
#			cls.PHE  = interp1d(fin["PHE" + "/r"], fin["PHE" + "/conductivity"], kind = "cubic")
#			cls.GLY  = interp1d(fin["GLY" + "/r"], fin["GLY" + "/conductivity"], kind = "cubic")
#			cls.HIS  = interp1d(fin["HIS" + "/r"], fin["HIS" + "/conductivity"], kind = "cubic")
#			cls.ILE  = interp1d(fin["ILE" + "/r"], fin["ILE" + "/conductivity"], kind = "cubic")
#			cls.LYS  = interp1d(fin["LYS" + "/r"], fin["LYS" + "/conductivity"], kind = "cubic")
#			cls.LEU  = interp1d(fin["LEU" + "/r"], fin["LEU" + "/conductivity"], kind = "cubic")
#			cls.MET  = interp1d(fin["MET" + "/r"], fin["MET" + "/conductivity"], kind = "cubic")
#			cls.ASN  = interp1d(fin["ASN" + "/r"], fin["ASN" + "/conductivity"], kind = "cubic")
#			cls.PRO  = interp1d(fin["PRO" + "/r"], fin["PRO" + "/conductivity"], kind = "cubic")
#			cls.GLN  = interp1d(fin["GLN" + "/r"], fin["GLN" + "/conductivity"], kind = "cubic")
#			cls.ARG  = interp1d(fin["ARG" + "/r"], fin["ARG" + "/conductivity"], kind = "cubic")
#			cls.SER  = interp1d(fin["SER" + "/r"], fin["SER" + "/conductivity"], kind = "cubic")
#			cls.THR  = interp1d(fin["THR" + "/r"], fin["THR" + "/conductivity"], kind = "cubic")
#			cls.VAL  = interp1d(fin["VAL" + "/r"], fin["VAL" + "/conductivity"], kind = "cubic")
#			cls.TRP  = interp1d(fin["TRP" + "/r"], fin["TRP" + "/conductivity"], kind = "cubic")
#			cls.TYR  = interp1d(fin["TYR" + "/r"], fin["TYR" + "/conductivity"], kind = "cubic")
#			cls.BULK = bulk_conductivity



#Residue2Conductivity = {
#	"DA"   : ConductivityOf.DA,
#	"DT"   : ConductivityOf.DT,
#	"DG"   : ConductivityOf.DG,
#	"DC"   : ConductivityOf.DC,
#	"ALA"  : ConductivityOf.ALA,
#	"CYS"  : ConductivityOf.CYS,
#	"ASP"  : ConductivityOf.ASP,
#	"GLU"  : ConductivityOf.GLU,
#	"PHE"  : ConductivityOf.PHE,
#	"GLY"  : ConductivityOf.GLY,
#	"HIS"  : ConductivityOf.HIS,
#	"ILE"  : ConductivityOf.ILE,
#	"LYS"  : ConductivityOf.LYS,
#	"LEU"  : ConductivityOf.LEU,
#	"MET"  : ConductivityOf.MET,
#	"ASN"  : ConductivityOf.ASN,
#	"PRO"  : ConductivityOf.PRO,
#	"GLN"  : ConductivityOf.GLN,
#	"ARG"  : ConductivityOf.ARG,
#	"SER"  : ConductivityOf.SER,
#	"THR"  : ConductivityOf.THR,
#	"VAL"  : ConductivityOf.VAL,
#	"TRP"  : ConductivityOf.TRP,
#	"TYR"  : ConductivityOf.TYR,
#	"BULK" : ConductivityOf.BULK
#}

Residue2Conductivity = dict()

def set_conductivity_parameters():

	if SEMParams.log_flag:
		print("Reading conductivity parameters from {} ...".format(SEMParams.conductivity_path))
	
	__fin = h5py.File(SEMParams.conductivity_path, 'r')
	Residue2Conductivity["DA"  ] = interp1d(__fin["DA"  + "/r"], __fin["DA"  + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["DT"  ] = interp1d(__fin["DT"  + "/r"], __fin["DT"  + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["DG"  ] = interp1d(__fin["DG"  + "/r"], __fin["DG"  + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["DC"  ] = interp1d(__fin["DC"  + "/r"], __fin["DC"  + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["ALA" ] = interp1d(__fin["ALA" + "/r"], __fin["ALA" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["CYS" ] = interp1d(__fin["CYS" + "/r"], __fin["CYS" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["ASP" ] = interp1d(__fin["ASP" + "/r"], __fin["ASP" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["GLU" ] = interp1d(__fin["GLU" + "/r"], __fin["GLU" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["PHE" ] = interp1d(__fin["PHE" + "/r"], __fin["PHE" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["GLY" ] = interp1d(__fin["GLY" + "/r"], __fin["GLY" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["HIS" ] = interp1d(__fin["HIS" + "/r"], __fin["HIS" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["ILE" ] = interp1d(__fin["ILE" + "/r"], __fin["ILE" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["LYS" ] = interp1d(__fin["LYS" + "/r"], __fin["LYS" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["LEU" ] = interp1d(__fin["LEU" + "/r"], __fin["LEU" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["MET" ] = interp1d(__fin["MET" + "/r"], __fin["MET" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["ASN" ] = interp1d(__fin["ASN" + "/r"], __fin["ASN" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["PRO" ] = interp1d(__fin["PRO" + "/r"], __fin["PRO" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["GLN" ] = interp1d(__fin["GLN" + "/r"], __fin["GLN" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["ARG" ] = interp1d(__fin["ARG" + "/r"], __fin["ARG" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["SER" ] = interp1d(__fin["SER" + "/r"], __fin["SER" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["THR" ] = interp1d(__fin["THR" + "/r"], __fin["THR" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["VAL" ] = interp1d(__fin["VAL" + "/r"], __fin["VAL" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["TRP" ] = interp1d(__fin["TRP" + "/r"], __fin["TRP" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["TYR" ] = interp1d(__fin["TYR" + "/r"], __fin["TYR" + "/conductivity"], kind = "cubic", bounds_error = False, fill_value = SEMConst.BULK),
	Residue2Conductivity["BULK"] = bulk_conductivity
	__fin.close()
