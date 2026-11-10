"""Plot points for the R-work / R-free curve on the results page.

The window draws with a Tk canvas, so all this module does is turn the
per-round metric rows of ``collect_results`` into percentages and pick a
readable y axis.  Keeping it toolkit-free lets the tests check the numbers.
"""

from __future__ import annotations


def round_series(metrics):
    """Return one row per refinement round that has an R value.

    Each row keeps the stage label, the round number (x axis) and the two
    R values as percentages; missing values stay ``None`` so a series can
    stop early without inventing points.
    """
    rows = []
    for row in metrics or []:
        work = row.get('r_work')
        free = row.get('r_free')
        if work is None and free is None:
            continue
        number = row.get('round')
        rows.append(dict(
            index=len(rows) + 1,
            round=int(number) if number is not None else len(rows) + 1,
            label=row.get('stage') or '第 %d 轮' % (len(rows) + 1),
            r_work=None if work is None else 100.0 * float(work),
            r_free=None if free is None else 100.0 * float(free),
        ))
    return rows


def y_range(rows):
    """Return a (low, high) window with a little air above and below."""
    values = [value for row in rows for value in (row['r_work'], row['r_free'])
              if value is not None]
    if not values:
        return 0.0, 50.0
    low, high = min(values), max(values)
    span = max(high - low, 4.0)
    return max(0.0, low - 0.15 * span), high + 0.15 * span


def y_ticks(low, high, count=5):
    """Round tick values (1/2/5 steps) that fit between ``low`` and ``high``."""
    if high <= low:
        return [low]
    steps = (1, 2, 5, 10, 20, 25, 50, 100)
    step = steps[-1]
    for candidate in steps:
        if (high - low) / candidate <= count:
            step = candidate
            break
    start = int(low // step) * step
    if start < low:
        start += step
    ticks = []
    value = start
    while value <= high + 1e-9:
        ticks.append(float(value))
        value += step
    if not ticks:
        ticks = [float(low), float(high)]
    if ticks[-1] < high:
        ticks.append(ticks[-1] + step)   # close the gap above the top value
    return ticks
