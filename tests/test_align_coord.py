from gbkviz.align_coord import AlignCoord


def test_is_inverted():
    """test is inverted"""
    align_coord = AlignCoord(11, 100, 501, 600, 90, 100, 80.0, "ref", "query")
    assert align_coord.is_inverted is False

    align_coord = AlignCoord(100, 11, 501, 600, 90, 100, 80.0, "ref", "query")
    assert align_coord.is_inverted is True


def test_add_offset():
    """test add offset"""
    ref_start, ref_end = 10, 100
    query_start, query_end = 500, 600
    align_coord = AlignCoord(
        ref_start, ref_end, query_start, query_end, 90, 100, 80.0, "ref", "query"
    )
    ref_offset, query_offset = 100, 150
    offset_align_coord = align_coord.add_offset(ref_offset, query_offset)
    assert (
        offset_align_coord.ref_start == ref_start + ref_offset
        and offset_align_coord.ref_end == ref_end + ref_offset
        and offset_align_coord.query_start == query_start + query_offset
        and offset_align_coord.query_end == query_end + query_offset
    )


def test_filter():
    """test filter"""
    align_coords = [
        AlignCoord(1, 1, 1, 1, 100, 100, 80.0, "ref", "query"),
        AlignCoord(1, 1, 1, 1, 150, 250, 90.0, "ref", "query"),
        AlignCoord(1, 1, 1, 1, 300, 300, 60.0, "ref", "query"),
    ]
    # No setting
    assert len(AlignCoord.filter(align_coords)) == 3
    # Min Length setting
    assert len(AlignCoord.filter(align_coords, min_length=130)) == 2
    assert len(AlignCoord.filter(align_coords, min_length=200)) == 1
    assert len(AlignCoord.filter(align_coords, min_length=500)) == 0
    # Identity setting
    assert len(AlignCoord.filter(align_coords, min_identity=70)) == 2
    assert len(AlignCoord.filter(align_coords, min_identity=85)) == 1
    assert len(AlignCoord.filter(align_coords, min_identity=95)) == 0
    # Both setting
    assert len(AlignCoord.filter(align_coords, 200, 70)) == 0
