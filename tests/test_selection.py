from phasebuild.selection import independent_set


def test_global_selection_can_replace_a_greedy_conflicting_fragment():
    # Fragment 0 conflicts with 1 and 2, which coexist and have a larger total score.
    selected = independent_set([6.,4.,4.], [6,3,3], [(1<<1)|(1<<2),1,1], capacity=6)
    assert selected == [1,2]


def test_selection_obeys_capacity_and_conflicts():
    selected = independent_set([6.,4.,4.], [6,3,3], [(1<<1)|(1<<2),1,1], capacity=5)
    assert selected in [[1],[2]]
