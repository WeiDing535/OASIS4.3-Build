from phasebuild import charts


def test_round_series_converts_to_percentages_and_skips_empty_rows():
    metrics = [
        {'stage': '第 1 轮', 'round': 1, 'r_work': 0.4103, 'r_free': 0.4389},
        {'stage': '第 2 轮', 'round': 2, 'r_work': 0.3791, 'r_free': 0.4022},
        {'stage': '第 3 轮', 'round': 3, 'r_work': None, 'r_free': None},
    ]
    rows = charts.round_series(metrics)
    assert [row['round'] for row in rows] == [1, 2]
    assert round(rows[0]['r_work'], 4) == 41.03
    assert round(rows[1]['r_free'], 4) == 40.22
    assert rows[1]['index'] == 2


def test_round_series_keeps_a_missing_r_free_as_none():
    rows = charts.round_series([{'stage': '第 1 轮', 'round': 1, 'r_work': 0.4, 'r_free': None}])
    assert rows[0]['r_work'] == 40.0
    assert rows[0]['r_free'] is None


def test_round_series_falls_back_to_the_row_position():
    rows = charts.round_series([{'r_work': 0.3}, {'r_work': 0.28}])
    assert [row['round'] for row in rows] == [1, 2]
    assert [row['label'] for row in rows] == ['第 1 轮', '第 2 轮']


def test_y_range_adds_air_and_never_goes_negative():
    low, high = charts.y_range([{'r_work': 31.0, 'r_free': 38.0}])
    assert low >= 0.0
    assert low < 31.0 and high > 38.0
    assert charts.y_range([]) == (0.0, 50.0)


def test_y_ticks_land_on_round_values_inside_the_window():
    ticks = charts.y_ticks(29.0, 45.0)
    assert ticks == [30.0, 35.0, 40.0, 45.0]
    assert all(29.0 <= tick <= 45.0 for tick in ticks)


def test_y_ticks_keep_the_axis_readable_for_a_wide_span():
    ticks = charts.y_ticks(10.0, 90.0, count=5)
    assert len(ticks) <= 6
    assert ticks[0] >= 10.0
    assert ticks[-1] <= 110.0          # at most one step above the top value
