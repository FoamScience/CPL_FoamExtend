#!/usr/bin/env python2
import pytest
import numpy as np
from cplpy import run_test, prepare_config, parametrize_file, get_test_dir
import os
import sys
import subprocess
import re
import json

try:
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
except:
    print "Error: PyFoam package is required to run the tests"
    sys.exit(1)

# -----Forces TESTS-----

# EXPLANATION:

MD_FNAME = "dummyMD_forces.py"
MD_ARGS = MD_FNAME
MD_EXEC = "python2"
CFD_FNAME = "test_forces_case"
CFD_ARGS = "-parallel -case " + CFD_FNAME
CFD_EXEC = "CPLIcoFoam"
TEST_TEMPLATE_DIR = os.path.join(get_test_dir(), "templates")
TEST_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture()
def prepare_config_fix(tmpdir):
    prepare_config(tmpdir, TEST_DIR, [MD_FNAME, CFD_FNAME, "config.cpl"])


def name_with_step(fname, step):
    split_fname = re.split('(\W)', fname)
    dot_pos = split_fname.index(".")
    split_fname.insert(dot_pos , str(step))
    return "".join(split_fname)

def compare_forces(tol, md_fname="md_forces.dat",
                   openfoam_casedir="test_forces_case", steps=1, dt=2, mode="test"):
    #NOTE: dt hast to be an integer for this to work
    for step in xrange(1, steps+1):
        openfoam_dir = os.path.join(openfoam_casedir, str(step*dt))
        # Get viscosity conversion
        viscosity_conversion = json.load(open("config.cpl",'r'))["conversion-factors"]["dyn-viscosity"]
        # OpenFOAM data   
        s_xy = ParsedParameterFile(openfoam_dir + "/sigmaxy")["internalField"]
        # WARNING: This wont work if the value is uniform
        s_yz = ParsedParameterFile(openfoam_dir + "/sigmayz")["internalField"]
        s_xz = ParsedParameterFile(openfoam_dir + "/sigmaxz")["internalField"]
        s_xx = ParsedParameterFile(openfoam_dir + "/sigmaxx")["internalField"]
        s_yy = ParsedParameterFile(openfoam_dir + "/sigmayy")["internalField"]
        s_zz = ParsedParameterFile(openfoam_dir + "/sigmazz")["internalField"]
        cell_cx = ParsedParameterFile(openfoam_dir + "/ccx")["internalField"]
        cell_cy = ParsedParameterFile(openfoam_dir + "/ccy")["internalField"]
        cell_cz = ParsedParameterFile(openfoam_dir + "/ccz")["internalField"]

        openfoam_cells = {}
        cell_no = 0
        for cell_no in xrange(0, len(s_xy)): 
            cell_coord = (float(cell_cx[cell_no]), float(cell_cy[cell_no]), float(cell_cz[cell_no]))
            # Openfoam output cell centres with 6 decimal figures
            k = "{0:.5f}".format(cell_coord[0]), "{0:.5f}".format(cell_coord[1]),\
                "{0:.5f}".format(cell_coord[2])
            openfoam_cells[k] = np.array([float(s_xy[cell_no]), float(s_yy[cell_no]), 
                                          float(s_yz[cell_no])])

        # Line format of dummy md forces file -- > x y z sxy syy szy
        md_fname_step =  name_with_step(md_fname, step)
        with open(md_fname_step, "r") as cfd_file:
            cfd_lines = cfd_file.readlines()
        md_lines = [l[:-1].split(" ") for l in cfd_lines]
        md_cells = {}

        # Get density
        density = json.load(open("config.cpl",'r'))["initial-conditions"]["density"]
        density_conversion = json.load(open("config.cpl",'r'))["conversion-factors"]["density"]
        for l in md_lines:
            k = "{0:.5f}".format(float(l[0])), "{0:.5f}".format(float(l[1])), "{0:.5f}".format(float(l[2]))
            md_cells[k] = np.array([float(l[3]), float(l[4]), float(l[5])])
        for k in md_cells.keys():
            try:
                openfoam_cells[k] *= density * density_conversion * viscosity_conversion
                diff_forces = abs(md_cells[k] - openfoam_cells[k])
                if (np.any(diff_forces > tol)):
                    print "Cell %s value differs in md : %s and cfd: %s" % (str(k), str(density*md_cells[k]), str(openfoam_cells[k]))
                    if mode == "test":
                        assert False
                    else:
                        sys.exit()
            except KeyError:
                print "Cell not found: cell " + str(k)
                if mode == "test":
                    assert False
                else:
                    sys.exit()
        print "Step %d: " %step, "SUCCESS"





# -----FORCES TESTS-----

# EXPLANATION: See README-test located in this folder.

@pytest.mark.parametrize("cfdprocs, mdprocs, cells, err_msg", [
                         ((3, 3, 3), (3, 3, 3), (15, 15, 15), ""),
                         ((1, 2, 1), (3, 2, 2), (30, 36, 24), ""),
                         ((4, 3, 3), (4, 3, 3), (20, 15, 27), ""),
                         ((3, 3, 3), (3, 3, 3), (30, 15, 21), "")])
def test_forcesC2P(prepare_config_fix, cfdprocs, mdprocs, cells,  err_msg):
    MD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0}
    MD_PARAMS["npx"], MD_PARAMS["npy"], MD_PARAMS["npz"] = mdprocs

    CFD_PARAMS = {"lx": 300.0, "ly": 300.0, "lz": 300.0}
    CFD_PARAMS["npx"], CFD_PARAMS["npy"], CFD_PARAMS["npz"] = cfdprocs
    CFD_PARAMS["ncx"], CFD_PARAMS["ncy"], CFD_PARAMS["ncz"] = cells
    # Needed for decomposParDict
    CFD_PARAMS["nprocs"] = cfdprocs[0]*cfdprocs[1]*cfdprocs[2]

    CONFIG_PARAMS = {"cfd_bcx": 1, "cfd_bcy": 0, "cfd_bcz": 0,
                     "olap_xlo": 1, "olap_xhi": cells[0],
                     "olap_ylo": 1, "olap_yhi": 5,
                     "olap_zlo": 1, "olap_zhi": cells[2],
                     "cnst_xlo": 1, "cnst_xhi": cells[0],
                     "cnst_ylo": 3, "cnst_yhi": 5,
                     "cnst_zlo": 1, "cnst_zhi": cells[2],
                     "bndry_xlo": 1, "bndry_xhi": cells[0],
                     "bndry_ylo": 1, "bndry_yhi": 1,
                     "bndry_zlo": 1, "bndry_zhi": cells[2],
                     "tstep_ratio": 1, }

    # Parametrize OpenFOAM files
    mesh_file = os.path.join(CFD_FNAME+"/", "constant/polyMesh/blockMeshDict") 
    parametrize_file(mesh_file, mesh_file, CFD_PARAMS)
    control_dict_file = os.path.join(CFD_FNAME+"/", "system/decomposeParDict") 
    parametrize_file(control_dict_file, control_dict_file, CFD_PARAMS)

    try:
        subprocess.check_output(["blockMesh", "-case", CFD_FNAME])
    except:
        assert False

    try:
        subprocess.check_output(["decomposePar", "-case", CFD_FNAME])
    except:
        assert False

    correct = run_test(TEST_TEMPLATE_DIR, CONFIG_PARAMS, MD_EXEC, MD_FNAME, MD_ARGS,
                       CFD_EXEC, CFD_FNAME, CFD_ARGS, MD_PARAMS, CFD_PARAMS, err_msg, True)
    if correct:
        # Reconstruct the fields from processor directories.
        try:
            subprocess.check_output(["reconstructPar", "-case", CFD_FNAME])
        except:
            assert False
        # Reconstruct the fields from processor directories.
        try:
            subprocess.check_output(["writeCellCentres", "-case", CFD_FNAME])
        except:
            assert False
        # Calculate stress components using ...
        try:
            subprocess.check_output(["stressComponents", "-case", CFD_FNAME])
        except:
            assert False
        compare_forces(1e-6, steps=2, dt=2, mode="test")
