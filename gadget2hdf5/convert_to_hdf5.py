import g3read
import h5py
import numpy as np

ptypes = [0, 1, 4, 5]


def convert_to_hdf5(infile, pos_name="POS ", vel_name="VEL "):

    infile = "uid0_r3_large_snap_v2box_031"
    outfile = f"{infile}.hdf5"

    f = g3read.GadgetFile(infile)

    field_map = {
        pos_name: "Coordinates",
        vel_name: "Velocities",
        "MASS": "Masses",
        "RHO ": "Density",
        "U   ": "InternalEnergy",
        "POT ": "Potential",
        "SFR ": "StarFormationRate",
        "HSML": "SmoothingLength",
        "TEMP": "Temperature",
        "HOTT": "HotTemperature",
        "AGE ": "StellarAge",
        "CLDX": "ColdFraction",
    }

    ff = h5py.File(outfile, "w")

    h = ff.create_group("Header")
    h.attrs["BoxSize"] = f.header.BoxSize
    h.attrs["Flag_Cooling"] = f.header.flag_cooling
    h.attrs["Flag_DoublePrecision"] = f.header.flag_doubleprecision
    h.attrs["Flag_Feedback"] = f.header.flag_feedback
    h.attrs["Flag_IC_Info"] = 3
    h.attrs["Flag_Sfr"] = f.header.flag_sfr
    h.attrs["Flag_StellarAge"] = f.header.flag_stellarage
    h.attrs["HubbleParam"] = f.header.HubbleParam
    h.attrs["MassTable"] = f.header.mass
    h.attrs["MassTable"][1] = 0.0
    h.attrs["NumFilesPerSnapshot"] = f.header.num_files
    h.attrs["NumPart_ThisFile"] = f.header.npart
    h.attrs["NumPart_Total"] = f.header.npartTotal
    h.attrs["NumPart_Total_HighWord"] = np.array([0] * 6, dtype="int32")
    h.attrs["Omega0"] = f.header.Omega0
    h.attrs["OmegaLambda"] = f.header.OmegaLambda
    h.attrs["Redshift"] = f.header.redshift
    h.attrs["Time"] = f.header.time

    elems = {
        11: ["He", "C", "Ca", "O", "N", "Ne", "Mg", "S", "Si", "Fe", "Ej"],
        15: [
            "He",
            "C",
            "Ca",
            "O",
            "N",
            "Ne",
            "Mg",
            "S",
            "Si",
            "Fe",
            "Na",
            "Al",
            "Ar",
            "Ni",
            "Ej",
        ],
    }

    gs = {i: ff.create_group(f"PartType{i}") for i in ptypes}

    for name, blk in f.blocks.items():
        for ptype in ptypes:
            if not blk.ptypes[ptype]:
                continue
            arr = f.read_new(name, ptype)
            g = gs[i]
            if name == "Zs  ":
                elem_size = arr.shape[1]
                if "Flag_Metals" not in h.attrs:
                    h.attrs["Flag_Metals"] = elem_size
                frac = arr / g["Masses"][()]
                g.create_dataset("Metallicity", data=frac[:, 1:].sum(axis=1))
                if "H" not in elems[elem_size]:
                    g.create_dataset("H_fraction", data=1.0 - frac.sum(axis=1))
                for i, elem in enumerate(elems[elem_size]):
                    g.create_dataset(f"{elem}_fraction", data=frac[:, i])
            elif name == "ID  ":
                g.create_dataset("ParticleIDs", data=arr)
            elif name in field_map:
                g.create_dataset(field_map[name], data=arr)

    ff.flush()
    ff.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert a Gadget Format 2 file to HDF5."
    )
    parser.add_argument("infile", type=str, help="The name of the file to convert.")
    parser.add_argument(
        "--pos_name",
        type=str,
        default="POS ",
        help="The name of the position field. Default: 'POS '.",
    )
    parser.add_argument(
        "--vel_name",
        type=str,
        default="VEL ",
        help="The name of the velocity field. Default: 'VEL '.",
    )
    args = parser.parse_args()

    convert_to_hdf5(args.infile, pos_name=args.pos_name, vel_name=args.vel_name)
