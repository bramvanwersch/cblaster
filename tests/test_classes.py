#!/usr/bin/env python3

"""
Test suite for classes.
"""

import pytest
from cblaster import classes


# HIT TESTS

@pytest.fixture
def hits():
    hits = [
        classes.Hit("q1", "s1", "70.90", "54.60", "0.0", "500.30"),
        classes.Hit("q1", "s2", "5.90", "21.60", "1.06e-29", "100.30"),
        classes.Hit("q2", "s3", "2.90", "34.60", "1.02e-30", "200.30"),
        classes.Hit("q3", "s5", "1123.90", "12.60", "0.01", "500.30"),
        classes.Hit("q4", "s1", "0.90", "1.60", "0.005", "500.30"),
        classes.Hit("q1", "s1", "70.90", "54.60", "0.0", "500.30")
    ]

    return hits


@pytest.fixture
def subjects(hits):
    subjects = [
        classes.Subject(hits[:2], "subject1", "ipg1", "0", "100", "-"),
        classes.Subject([hits[1], hits[0]], "subject1", "ipg1", 0, 100, "-"),
        classes.Subject(hits[2:3], "subject2", "ipg2", 1200, 6000, "+"),
        classes.Subject(hits[3:4], "subject3", "ipg3", 500, "1000", "+"),
        classes.Subject(hits[4:5], "subject4", "ipg4", 6453, 8000, "+"),
    ]
    return subjects


@pytest.fixture
def clusters(subjects):
    clusters = [
        classes.Cluster([0, 1], subjects[:2]),
        classes.Cluster([0, 1], subjects[:2], 100, 50, 200),
        classes.Cluster([0, 2, 3, 4], [subjects[0], *subjects[2:]]),
        classes.Cluster([0, 1], [subjects[1], subjects[0]])
    ]
    return clusters


@pytest.fixture
def scaffolds(subjects, clusters):
    scaffolds = [
        classes.Scaffold("scaffold1", [clusters[2], clusters[0]], subjects[:5]),
        classes.Scaffold("scaffold2", [clusters[0]], subjects[:2])
    ]

    return scaffolds


@pytest.fixture
def scaffold_summaries():
    with open("data/scaffold_summaries") as f:
        summary_text = f.read()
    return [s for s in summary_text.split("\n>>\n") if len(s) > 0]


@pytest.fixture
def organisms(scaffolds):
    organisms = [
        classes.Organism("organism1", "strain_organism1", {scaffold.accession: scaffold for scaffold in scaffolds[0:1]})
    ]
    return organisms


@pytest.fixture
def organism_summaries():
    with open("data/organism_summaries") as f:
        summary_text = f.read()
    return [s for s in summary_text.split("\n>>\n") if len(s) > 0]


def test_hit_instantiation(hits):
    assert isinstance(hits[0], classes.Serializer)
    assert [str, str, float, float, float, float] == [
        type(getattr(hits[0], val))
        for val in [
            "query",
            "subject",
            "identity",
            "coverage",
            "evalue",
            "bitscore"
        ]
    ]


@pytest.mark.parametrize(
    "index, decimals, result",
    [
        (0, 4, ["q1", "s1", "70.9", "54.6", "0", "500.3"]),
        (0, 0, ["q1", "s1", "71", "55", "0", "500"]),
        (1, 2, ["q1", "s2", "5.9", "21.6", "1.1e-29", "100.3"]),
    ],
)
def test_hit_values(hits, index, decimals, result):
    assert hits[index].values(decimals) == result


def test_hit_str(hits):
    # this seems to be incorrect at the moment. The coverage and identity are formatted
    # into a percentage, when they are already percentage values
    assert str(hits[0]) == "Hit: q1 - s1: 7090.00%/5460.00%"


def test_hit_equals(hits):
    with pytest.raises(NotImplementedError):
        bool_ = hits[0] == "invalid value"
    assert hits[0] == hits[5]
    assert hits[0] != hits[4]


def test_hit_copy(hits):
    copy_hit = hits[0].copy(test="test_value")
    assert ["q1", "s1", 70.9, 54.6, 0.0, 500.3, "test_value"] == [
        getattr(copy_hit, val) for val in [
            "query",
            "subject",
            "identity",
            "coverage",
            "evalue",
            "bitscore",
            "test"]
    ]


def test_hit_to_dict(hits):
    assert hits[0].to_dict() == {
        "query": "q1",
        "subject": "s1",
        "identity": 70.9,
        "coverage": 54.6,
        "evalue": 0.0,
        "bitscore": 500.3
    }


def test_hit_from_dict():
    from_dict_class = classes.Hit.from_dict({
        "query": "q1",
        "subject": "s1",
        "identity": 70.9,
        "coverage": 54.6,
        "evalue": 0.0,
        "bitscore": 500.3})
    assert ["q1", "s1", 70.9, 54.6, 0.0, 500.3] == [
        getattr(from_dict_class, val) for val in [
            "query",
            "subject",
            "identity",
            "coverage",
            "evalue",
            "bitscore",
        ]
    ]


# SUBJECT TESTS

def test_subject_instantiation(subjects):
    assert isinstance(subjects[0], classes.Serializer)
    assert [list, str, str, int, int, str] == [
        type(getattr(subjects[0], val))
        for val in [
            "hits",
            "ipg",
            "name",
            "start",
            "end",
            "strand"
        ]
    ]


def test_subject_empty_instantiation():
    assert [[], None, None, None, None, None] == [
        getattr(classes.Subject(), val)
        for val in [
            "hits",
            "ipg",
            "name",
            "start",
            "end",
            "strand"
        ]
    ]


def test_subject_equals(subjects):
    with pytest.raises(NotImplementedError):
        bool_ = subjects[0] == "invalid value"
    assert subjects[0] == subjects[1]
    assert subjects[0] != subjects[2]


def test_subject_to_dict(hits, subjects):
    assert subjects[3].to_dict() == {
        "hits": [hits[3].to_dict()],
        "name": "subject3",
        "ipg": "ipg3",
        "start": 500,
        "end": 1000,
        "strand": "+",
    }


@pytest.mark.parametrize(
    "index, decimals, result",
    [
        (0, 4, [("q1", "s1", "70.9", "54.6", "0", "500.3", "0", "100", "-"),
                ("q1", "s2", "5.9", "21.6", "1.06e-29", "100.3", "0", "100", "-")]),
        (0, 0, [("q1", "s1", "71", "55", "0", "500", "0", "100", "-"),
                ("q1", "s2", "6", "22", "1e-29", "100", "0", "100", "-")])
    ],
)
def test_subject_values(subjects, index, decimals, result):
    assert subjects[index].values(decimals) == result


def test_subject_from_dict(hits, subjects):
    from_dict_class = classes.Subject.from_dict(
        {
            "hits": [hits[3].to_dict()],
            "name": "subject3",
            "ipg": "ipg3",
            "start": 500,
            "end": 1000,
            "strand": "+",
        })
    assert [[hits[3]], "subject3", "ipg3", 500, 1000, "+"] == [
        getattr(from_dict_class, val) for val in [
            "hits",
            "name",
            "ipg",
            "start",
            "end",
            "strand",
        ]
    ]


# CLUSTER TESTS


def test_cluster_instantiation(clusters):
    assert isinstance(clusters[0], classes.Serializer)
    assert [list, list, float, int, int] == [
        type(getattr(clusters[0], val)) for val in [
            "indices",
            "subjects",
            "score",
            "start",
            "end",
        ]
    ]


def test_cluster_empty_instantiation():
    assert [[], [], 0, None, None] == [
        getattr(classes.Cluster(), val) for val in [
            "indices",
            "subjects",
            "score",
            "start",
            "end",
        ]
    ]


def test_cluster_iter(subjects, clusters):
    assert subjects[:2] == [s for s in clusters[0]]


def test_cluster_len(clusters):
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2
    assert len(clusters[2]) == 4


def test_cluster_equals(clusters):
    with pytest.raises(NotImplementedError):
        bool_ = clusters[0] == "invalid value"
    assert clusters[0] != clusters[1]
    assert clusters[0] == clusters[3]


def test_cluster_calculate_score(clusters):
    assert clusters[0].calculate_score() == 2.10006
    assert clusters[2].calculate_score() == 4.17012
    assert clusters[2].calculate_score(["q1", "q3", "q2", "q4"]) == 5.17012
    assert clusters[2].calculate_score(["q1", "q2", "q3", "q4"]) == 7.17012


def test_cluster_to_dict(clusters):
    assert clusters[0].to_dict() == {
        "indices": [0, 1],
        "score": 2.10006,
        "start": 0,
        "end": 100,
    }


def test_cluster_from_dict(subjects, clusters):
    from_dict_class = classes.Cluster.from_dict({
        "indices": [0, 1],
        "score": 2.10006,
        "start": 0,
        "end": 100,
    },
        *subjects[:2]
    )
    assert [[0, 1], subjects[:2], 2.10006, 0, 100] == [
        getattr(from_dict_class, val) for val in [
            "indices",
            "subjects",
            "score",
            "start",
            "end"
        ]
    ]


# TEST SCAFFOLD

def test_scaffold_instantiation(scaffolds):
    assert isinstance(scaffolds[0], classes.Serializer)
    assert [str, list, list] == [
        type(getattr(scaffolds[0], val)) for val in [
            "accession",
            "subjects",
            "clusters"
        ]
    ]


def test_scaffold_empty_instantiation():
    assert ["test_name", [], []] == [
        getattr(classes.Scaffold("test_name"), val) for val in [
            "accession",
            "subjects",
            "clusters"
        ]
    ]


def test_scaffold_string(scaffolds):
    assert str(scaffolds[0]) == "SCAFFOLD: scaffold1 [5 hits in 2 clusters]"


def test_scaffold_add_clusters(subjects, scaffolds):
    scaffolds[0].add_clusters([subjects[1:4]])
    assert len(scaffolds[0].clusters) == 3
    assert all(subject in subjects[1:4] for subject in scaffolds[0].clusters[1].subjects)


@pytest.mark.parametrize(
    "index, hide_headers, delimiter, decimals, result_index",
    [
        (0, False, None, 4, 0),
        (0, True, None, 4, 1),
        (0, False, "____", 4, 2),
        (0, False, None, 0, 3)
    ],
)
def test_scaffold_summary(scaffolds, index, hide_headers, delimiter, decimals, result_index, scaffold_summaries):
    assert scaffolds[index].summary(hide_headers, delimiter, decimals) == scaffold_summaries[result_index]


def test_scaffold_to_dict(scaffolds):
    assert scaffolds[1].to_dict() == {
        'accession': 'scaffold2',
        'clusters': [cluster.to_dict() for cluster in scaffolds[1].clusters],
        'subjects': [subject.to_dict() for subject in scaffolds[1].subjects]
    }


def test_scaffold_from_dict(subjects, clusters, scaffolds):
    from_dict_class = classes.Scaffold.from_dict({
        'accession': 'scaffold2',
        'clusters': [cluster.to_dict() for cluster in scaffolds[1].clusters],
        'subjects': [subject.to_dict() for subject in scaffolds[1].subjects]
    })
    assert ["scaffold2", [clusters[0]], subjects[:2]] == [
        getattr(from_dict_class, val) for val in [
            "accession",
            "clusters",
            "subjects"
        ]
    ]


# TEST ORGANISMS

def test_organism_instantiation(organisms):
    assert isinstance(organisms[0], classes.Serializer)
    assert [str, str, dict] == [
        type(getattr(organisms[0], val)) for val in [
            "name",
            "strain",
            "scaffolds"
        ]
    ]


def test_organism_empty_instantiation():
    assert ["test_name", "test_strain", {}] == [
        getattr(classes.Organism("test_name", "test_strain"), val) for val in [
            "name",
            "strain",
            "scaffolds"
        ]
    ]


def test_organism_string(organisms):
    assert str(organisms[0]) == "ORGANISM: organism1 strain_organism1 [5 subjects on 1 scaffolds]"


def test_organism_clusters(clusters, organisms):
    assert organisms[0].clusters == [clusters[2], clusters[0]]


def test_organism_total_hit_clusters(organisms):
    assert organisms[0].total_hit_clusters == 2


@pytest.mark.parametrize(
    "index, hide_headers, delimiter, decimals, result_index",
    [
        (0, False, None, 4, 0),
        (0, True, None, 4, 1),
        (0, False, "____", 4, 2),
        (0, False, None, 0, 3)
    ],
)
def test_organism_summary(organisms, index, hide_headers, delimiter, decimals, result_index, organism_summaries):
    assert organisms[index].summary(hide_headers=hide_headers, delimiter=delimiter, decimals=decimals) ==\
           organism_summaries[result_index]


def test_organism_full_name(organisms):
    pass