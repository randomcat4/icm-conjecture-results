#!/usr/bin/env python3
"""Compare the candidate's serialized claims to independently generated data."""

from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)
CANDIDATE = HERE / "exact_rational_1e4.json"


def frac(s):
    return F(s)


def main():
    independent = json.loads((HERE / "independent_results.json").read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))

    assert candidate["dimension"] == 6
    assert candidate["common_denominator"] == 10000
    assert candidate["K_integer_numerators"] == [
        [5000, 0, 0, 3298, -1672, 3298],
        [0, 5000, 0, -1672, 3298, 3298],
        [0, 0, 5000, 3298, 3298, -1672],
        [3298, -1672, 3298, 5000, 0, 0],
        [-1672, 3298, 3298, 0, 5000, 0],
        [3298, 3298, -1672, 0, 0, 5000],
    ]
    assert candidate["B_integer_numerators"] == [
        [0, 11, -11, 0, 0, 0],
        [-11, 0, 11, 0, 0, 0],
        [11, -11, 0, 0, 0, 0],
        [0, 0, 0, 0, -11, 11],
        [0, 0, 0, 11, 0, -11],
        [0, 0, 0, -11, 11, 0],
    ]

    assert len(candidate["class_table"]) == len(independent["classes"]) == 8
    for got, expected in zip(candidate["class_table"], independent["classes"]):
        assert got["multiplicity"] == expected["multiplicity"]
        assert got["all_masks_decimal"] == expected["masks"]
        assert frac(got["center_probability"]) == frac(expected["center_probability"])
        assert frac(got["endpoint_probability"]) == frac(expected["endpoint_probability"])

    assert frac(candidate["minimum_center_probability"]) == frac(independent["minimum_atoms"]["center"])
    assert frac(candidate["minimum_endpoint_probability"]) == frac(independent["minimum_atoms"]["plus"])
    assert candidate["center_probabilities_sum"] == "1"
    assert candidate["endpoint_probabilities_sum"] == "1"
    assert candidate["endpoint_laws_equal_eventwise"] is True

    candidate_gap = (frac(candidate["gap_lower_fraction"]), frac(candidate["gap_upper_fraction"]))
    independent_gap = tuple(frac(x) for x in independent["gap_interval"])
    assert candidate_gap[0] <= independent_gap[0] <= independent_gap[1] <= candidate_gap[1]
    assert candidate_gap[0] > F(1, 3_700_000)
    assert candidate_gap[1] < F(1, 3_600_000)
    assert candidate_gap[1] < F(1, 1_000_000)

    print(json.dumps({
        "all_eight_candidate_classes_match_independent_output": True,
        "all_candidate_class_masks_match": True,
        "candidate_gap_contains_independent_tighter_gap": True,
        "candidate_simple_gap_bounds_valid": True,
        "frozen_primary_threshold_met": False,
        "status": "PASS",
    }, indent=2))


if __name__ == "__main__":
    main()
