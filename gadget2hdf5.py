import g3read
import h5py
import numpy as np

infile = "uid0_r3_large_snap_v2box_031"
outfile = f"{infile}.hdf5"

f = g3read.GadgetFile(infile)

field_map = {"XPOS": "Coordinates",
             "VPOS": "Velocities",
             "MASS": "Masses",
             "RHO ": "Density",
             "U   ": "InternalEnergy",
             "POT ": "Potential",
             "SFR ": "StarFormationRate",
             "HSML": "SmoothingLength",
             "TEMP": "Temperature",
             "HOTT": "HotTemperature",
             "AGE ": "StellarAge",
             "CLDX": "ColdFraction"}

ff = h5py.File(outfile, "w")

h = ff.create_group("Header")
h.attrs["BoxSize"] = f.header.BoxSize
h.attrs["Flag_Cooling"] = f.header.flag_cooling
h.attrs["Flag_DoublePrecision"] = f.header.flag_doubleprecision
h.attrs["Flag_Feedback"] = f.header.flag_feedback
h.attrs["Flag_IC_Info"] = 3
h.attrs["Flag_Metals"] = -1 # Will be set later
h.attrs["Flag_Sfr"] = f.header.flag_sfr
h.attrs["Flag_StellarAge"] = f.header.flag_stellarage
h.attrs["HubbleParam"] = f.header.HubbleParam
h.attrs["MassTable"] = f.header.mass
h.attrs["MassTable"][1] = 0.0
h.attrs["NumFilesPerSnapshot"] = f.header.num_files
h.attrs["NumPart_ThisFile"] = f.header.npart
h.attrs["NumPart_Total"] = f.header.npartTotal
h.attrs["NumPart_ThisFile"][1:] = 0
h.attrs["NumPart_Total"][1:] = 0
h.attrs["NumPart_Total_HighWord"] = np.array([0]*6, dtype='int32')
h.attrs["Omega0"] = f.header.Omega0
h.attrs["OmegaLambda"] = f.header.OmegaLambda
h.attrs["Redshift"] = f.header.redshift
h.attrs["Time"] = f.header.time

elems = {
    11: ["He", "C", "Ca", "O", "N", "Ne", "Mg", "S", "Si", "Fe", "Ej"],
    15: ["He", "C", "Ca", "O", "N", "Ne", "Mg", 
         "S", "Si", "Fe", "Na", "Al", "Ar", "Ni", "Ej"]
}

ggas = ff.create_group("PartType0")

for name, blk in f.blocks.items():
    if not blk.ptypes[0]:
        continue
    arr = f.read_new(name, 0)
    if name == "Zs  ":
        elem_size = arr.shape[1]
        h.attrs["Flag_Metals"] = elem_size
        mass = ggas["Masses"][()]
        ggas.create_dataset("Metallicity", data=arr[:,1:].sum(axis=1)/mass)
        if "H" not in elems[elem_size]:
            ggas.create_dataset("H_fraction", data=(1.0-arr.sum(axis=1)/mass))
        for i, elem in enumerate(elems[elem_size]):
            ggas.create_dataset(f"{elem}_fraction", data=arr[:,i]/mass)
    elif name == "ID  ":
        ggas.create_dataset("particle_index", data=arr)
    elif name in field_map:
        ggas.create_dataset(field_map[name], data=arr)

gdm = ff.create_group("PartType1")

for name, blk in f.blocks.items():
    if not blk.ptypes[1]:
        continue
    arr = f.read_new(name, 1)
    if name == "ID  ":
        gdm.create_dataset("particle_index", data=arr)
    elif name in field_map:
        gdm.create_dataset(field_map[name], data=arr)

gs = ff.create_group("PartType4")

for name, blk in f.blocks.items():
    if not blk.ptypes[4]:
        continue
    arr = f.read_new(name, 4)
    if name == "Zs  ":
        elem_size = arr.shape[1]
        mass = gs["Masses"][()]
        gs.create_dataset("Metallicity", data=arr[:,1:].sum(axis=1)/mass)
        if "H" not in elems[elem_size]:
            gs.create_dataset("H_fraction", data=(1.0-arr.sum(axis=1)/mass))
        for i, elem in enumerate(elems[elem_size]):
            gs.create_dataset(f"{elem}_fraction", data=arr[:,i]/mass)
    elif name == "ID  ":
        gs.create_dataset("particle_index", data=arr)
    elif name in field_map:
        gs.create_dataset(field_map[name], data=arr)
        
ff.close()

