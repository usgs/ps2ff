
import pandas as pd
import shutil
from impactutils.io.cmd import get_command_output


def test_single():
    cmd = "bin/run_ps2ff_single_event tests/config/test_single.ini"
    rc, so, se = get_command_output(cmd)

    rjb1 = pd.DataFrame.from_csv(
            "tests/data/Rjb_bytheta_Ratios.csv",
            header=6)
    vjb1 = pd.DataFrame.from_csv(
            "tests/data/Rjb_bytheta_Var.csv",
            header=6)
    rjb2 = pd.DataFrame.from_csv(
            "DataSingle/Rjb_bytheta_Ratios.csv",
            header=6)
    vjb2 = pd.DataFrame.from_csv(
            "DataSingle/Rjb_bytheta_Var.csv",
            header=6)

    rrup1 = pd.DataFrame.from_csv(
            "tests/data/Rrup_bytheta_Ratios.csv",
            header=6)
    vrup1 = pd.DataFrame.from_csv(
            "tests/data/Rrup_bytheta_Var.csv",
            header=6)
    rrup2 = pd.DataFrame.from_csv(
            "DataSingle/Rrup_bytheta_Ratios.csv",
            header=6)
    vrup2 = pd.DataFrame.from_csv(
            "DataSingle/Rrup_bytheta_Var.csv",
            header=6)

    pd.util.testing.assert_frame_equal(rjb1, rjb2)
    pd.util.testing.assert_frame_equal(vjb1, vjb2)
    pd.util.testing.assert_frame_equal(rrup1, rrup2)
    pd.util.testing.assert_frame_equal(vrup1, vrup2)

    # Clean up
    shutil.rmtree('DataSingle')


def test_single_by_theta_false():
    cmd = "bin/run_ps2ff_single_event tests/config/test_single_by_theta_false.ini"
    rc, so, se = get_command_output(cmd)

    rjb1 = pd.DataFrame.from_csv(
            "tests/data/Rjb_Ratios_by_theta_false.csv",
            header=6)
    vjb1 = pd.DataFrame.from_csv(
            "tests/data/Rjb_Var_by_theta_false.csv",
            header=6)
    rjb2 = pd.DataFrame.from_csv(
            "DataSingle/Rjb_Ratios.csv",
            header=6)
    vjb2 = pd.DataFrame.from_csv(
            "DataSingle/Rjb_Var.csv",
            header=6)

    rrup1 = pd.DataFrame.from_csv(
            "tests/data/Rrup_Ratios_by_theta_false.csv",
            header=6)
    vrup1 = pd.DataFrame.from_csv(
            "tests/data/Rrup_Var_by_theta_false.csv",
            header=6)
    rrup2 = pd.DataFrame.from_csv(
            "DataSingle/Rrup_Ratios.csv",
            header=6)
    vrup2 = pd.DataFrame.from_csv(
            "DataSingle/Rrup_Var.csv",
            header=6)

    pd.util.testing.assert_frame_equal(rjb1, rjb2)
    pd.util.testing.assert_frame_equal(vjb1, vjb2)
    pd.util.testing.assert_frame_equal(rrup1, rrup2)
    pd.util.testing.assert_frame_equal(vrup1, vrup2)

    # Clean up
    shutil.rmtree('DataSingle')


def test_single_N1():
    cmd = "bin/run_ps2ff_single_event tests/config/test_single_N1.ini"
    rc, so, se = get_command_output(cmd)

    rjb1 = pd.DataFrame.from_csv(
            "tests/data/Rjb_bytheta_Ratios_N1.csv",
            header=6)
    vjb1 = pd.DataFrame.from_csv(
            "tests/data/Rjb_bytheta_Var_N1.csv",
            header=6)
    rjb2 = pd.DataFrame.from_csv(
            "DataSingle/Rjb_bytheta_Ratios.csv",
            header=6)
    vjb2 = pd.DataFrame.from_csv(
            "DataSingle/Rjb_bytheta_Var.csv",
            header=6)

    rrup1 = pd.DataFrame.from_csv(
            "tests/data/Rrup_bytheta_Ratios_N1.csv",
            header=6)
    vrup1 = pd.DataFrame.from_csv(
            "tests/data/Rrup_bytheta_Var_N1.csv",
            header=6)
    rrup2 = pd.DataFrame.from_csv(
            "DataSingle/Rrup_bytheta_Ratios.csv",
            header=6)
    vrup2 = pd.DataFrame.from_csv(
            "DataSingle/Rrup_bytheta_Var.csv",
            header=6)

    pd.util.testing.assert_frame_equal(rjb1, rjb2)
    pd.util.testing.assert_frame_equal(vjb1, vjb2)
    pd.util.testing.assert_frame_equal(rrup1, rrup2)
    pd.util.testing.assert_frame_equal(vrup1, vrup2)

    # Clean up
    shutil.rmtree('DataSingle')
