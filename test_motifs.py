'''

'''
from nose.tools import nottest

from motifs import count_outgoing, count_incoming, all_nodes_where_incoming_sum_is_x_are_mech_y, gen_unique_cms
from generate_tpm import OR, AND

def test_count_outgoing():
    cm = [[0,1,1],
          [1,0,0],
          [0,0,0]]
    assert count_outgoing(cm) == [2,1,0]

def test_count_incoming():
    cm = [[0,0,1],
          [1,0,0],
          [1,0,0]]
    assert count_incoming(cm) == [2,0,1]


def test_all_nodes_where_incoming_sum_is_x_are_mech_y():
    cm = [[0,0,1],
          [1,0,0],
          [1,0,0]]
    nodes1 = [OR, OR, AND]
    nodes2 = [OR, OR, OR]

    assert not all_nodes_where_incoming_sum_is_x_are_mech_y(cm, nodes1, 1, [OR])
    assert all_nodes_where_incoming_sum_is_x_are_mech_y(cm, nodes2, 1, [OR])


@nottest
def test_gen_unique_cms():
    # A003085 - http://oeis.org/A003085
    NUM_NODES = 3
    expected = []

    actual = gen_unique_cms(NUM_NODES)

    assert expected == actual
