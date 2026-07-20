from app import branding


def test_slice_meter_html_marks_completed_slices():
    html = branding.slice_meter_html(2, 5, caption="Setup progress 2/5")
    # one wrapper per total slice
    assert html.count("mango-slice-wrap") == 5
    # exactly two filled (ripe) slices use the ripe flesh color
    assert html.count(branding._FLESH) == 2
    # the most recently filled slice gets the pop animation class
    assert html.count("mango-slice-wrap new") == 1
    assert "Setup progress 2/5" in html


def test_slice_meter_html_clamps_bounds():
    none_done = branding.slice_meter_html(0, 5)
    assert none_done.count("mango-slice-wrap") == 5
    assert none_done.count("mango-slice-wrap new") == 0

    over = branding.slice_meter_html(9, 5)
    assert over.count("mango-slice-wrap") == 5
    assert over.count(branding._FLESH) == 5  # completed clamped to total

    negative = branding.slice_meter_html(-3, 4)
    assert negative.count(branding._FLESH) == 0


def test_slicing_loader_contains_animation_and_label():
    html = branding.slicing_loader("Working…")
    assert "mango-loader" in html
    assert "mango-slice" in html
    assert "Working…" in html


def test_accent_html_renders_ripe_slices():
    assert branding.accent_html(3).count("mango-slice") == 3
    assert branding.accent_html(0).count("mango-slice") == 0


def test_style_defines_keyframes():
    assert "@keyframes mango-rise" in branding._STYLE
    assert "@keyframes mango-pop" in branding._STYLE
    assert "@keyframes mango-slicing" in branding._STYLE
