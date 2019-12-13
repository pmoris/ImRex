from src.data.controlCDR3Source import ControlCDR3Source


def test_control_cdr3():
    cdr3_control = ControlCDR3Source()

    # assert size matches de-duplicated size
    # NOTE: might differ from values seen elsewhere, since a variable size restriction
    #       can still be applied in other places...
    assert cdr3_control.data.shape == (487_437, 3)

    # assert size matches de-duplicated size after standard 10:20 size filter
    assert cdr3_control.data.loc[
        (cdr3_control.data["CDR3_beta"].str.len() >= 10)
        & (cdr3_control.data["CDR3_beta"].str.len() <= 20)
    ].shape == (484_914, 3)

    assert cdr3_control.data.loc[cdr3_control.data.duplicated("CDR3_beta")].size == 0
