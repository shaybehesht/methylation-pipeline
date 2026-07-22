import pysam

from core.qc import any_mn_tag, has_mn_tag


def _write_bam(path, *, with_mn: bool) -> str:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        read = pysam.AlignedSegment(bam.header)
        read.query_name = "r1"
        read.query_sequence = "ACGTACGT"
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 10
        read.mapping_quality = 60
        read.cigartuples = [(0, 8)]
        if with_mn:
            read.set_tag("MN", 8)
        bam.write(read)
    return str(path)


def test_has_mn_tag_true(tmp_path):
    bam = _write_bam(tmp_path / "mn.bam", with_mn=True)
    assert has_mn_tag(bam) is True


def test_has_mn_tag_false(tmp_path):
    bam = _write_bam(tmp_path / "no_mn.bam", with_mn=False)
    assert has_mn_tag(bam) is False


def test_has_mn_tag_unreadable(tmp_path):
    missing = tmp_path / "nope.bam"
    assert has_mn_tag(str(missing)) is None


def test_any_mn_tag_pacbio_like_all_missing(tmp_path):
    bams = [
        _write_bam(tmp_path / "a.bam", with_mn=False),
        _write_bam(tmp_path / "b.bam", with_mn=False),
        _write_bam(tmp_path / "c.bam", with_mn=False),
    ]
    # every BAM readable and none carry MN -> combine-strands should be disabled
    assert any_mn_tag(bams) is False


def test_any_mn_tag_ont_like_present(tmp_path):
    bams = [
        _write_bam(tmp_path / "a.bam", with_mn=False),
        _write_bam(tmp_path / "b.bam", with_mn=True),
    ]
    assert any_mn_tag(bams) is True


def test_any_mn_tag_unreadable_defaults_true(tmp_path):
    # can't determine -> don't override the user's choice
    assert any_mn_tag([str(tmp_path / "missing.bam")]) is True
