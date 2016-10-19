
import pandas as pd
import shutil
from impactutils.io.cmd import get_command_output


def test_rjb():
    cmd = "progs/RjbMeanVar.py progs/test_Rjb.ini"
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv("tests/data/test_Rjb_WC94_mechA_ar1p7_seis0_20_Ratios.csv", header = 6)
    v1 = pd.DataFrame.from_csv("tests/data/test_Rjb_WC94_mechA_ar1p7_seis0_20_Var.csv", header = 6)
    r2 = pd.DataFrame.from_csv("TestData/test_Rjb_WC94_mechA_ar1p7_seis0_20_Ratios.csv", header = 6)
    v2 = pd.DataFrame.from_csv("TestData/test_Rjb_WC94_mechA_ar1p7_seis0_20_Var.csv", header = 6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')

def test_rrup():
    cmd = "progs/RrupMeanVar.py progs/test_Rrup.ini"
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv("tests/data/test_Rrup_S14_mechA_ar2p0_seis0_15_Ratios.csv", header = 6)
    v1 = pd.DataFrame.from_csv("tests/data/test_Rrup_S14_mechA_ar2p0_seis0_15_Var.csv", header = 6)
    r2 = pd.DataFrame.from_csv("TestData/test_Rrup_S14_mechA_ar2p0_seis0_15_Ratios.csv", header = 6)
    v2 = pd.DataFrame.from_csv("TestData/test_Rrup_S14_mechA_ar2p0_seis0_15_Var.csv", header = 6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')

