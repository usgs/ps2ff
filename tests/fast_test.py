
import pandas as pd
import shutil
from impactutils.io.cmd import get_command_output


def test_rrup_HB08():
    conf = 'fast_rrup_HB08.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/Rrup_HB08_mechA_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/Rrup_HB08_mechA_ar1p0_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rrup_HB08_mechA_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rrup_HB08_mechA_ar1p0_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rrup_WC94_SS_F():
    conf = 'fast_rrup_WC94_SS_F.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechSS_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechSS_ar1p0_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechSS_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechSS_ar1p0_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rjb_WC94_SS_F():
    conf = 'fast_rjb_WC94_SS_F.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechSS_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechSS_ar1p0_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechSS_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechSS_ar1p0_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rrup_WC94_SS_T():
    conf = 'fast_rrup_WC94_SS_T.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechSS_LW_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechSS_LW_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechSS_LW_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechSS_LW_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rjb_WC94_SS_T():
    conf = 'fast_rjb_WC94_SS_T.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechSS_LW_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechSS_LW_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechSS_LW_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechSS_LW_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rrup_WC94_R_T():
    conf = 'fast_rrup_WC94_R_T.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechR_LW_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechR_LW_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechR_LW_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechR_LW_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rjb_WC94_R_T():
    conf = 'fast_rjb_WC94_R_T.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechR_LW_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechR_LW_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechR_LW_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechR_LW_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rrup_WC94_R_F():
    conf = 'fast_rrup_WC94_R_F.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechR_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rrup_WC94_mechR_ar1p0_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechR_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rrup_WC94_mechR_ar1p0_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')


def test_rjb_WC94_R_F():
    conf = 'fast_rjb_WC94_R_F.ini'
    cmd = "run_ps2ff tests/config/%s" % conf
    rc, so, se = get_command_output(cmd)
    r1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechR_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v1 = pd.DataFrame.from_csv(
        "tests/data/test_Rjb_WC94_mechR_ar1p0_seis0_15_Var.csv",
        header=6)
    r2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechR_ar1p0_seis0_15_Ratios.csv",
        header=6)
    v2 = pd.DataFrame.from_csv(
        "TestData/Rjb_WC94_mechR_ar1p0_seis0_15_Var.csv",
        header=6)

    pd.util.testing.assert_frame_equal(r1, r2)
    pd.util.testing.assert_frame_equal(v1, v2)

    # Clean up
    shutil.rmtree('TestData')

