from plot_embedding import get_nsc_id_from_path


def test_1():
    path = '/data/nci_open_training_data/191910.gexf.g2v3'
    assert 191910 == get_nsc_id_from_path(path)
