"""Spec 035 D2: AutoDock4 grid points use the even ceiling (AutoGrid needs an even npts)."""
import math

import pytest

from docking.runners.autodock4 import AutoDock4Runner


def test_d2_example_twenty_angstrom_box_at_point_375():
    # D2: 20 A at 0.375 A gives 54 points (round() gave 53, an odd count AutoGrid rejects).
    assert AutoDock4Runner._grid_points(20.0, 0.375) == 54


@pytest.mark.parametrize("size", [4.0, 9.7, 20.0, 21.0, 22.5, 31.1])
def test_grid_points_are_even_and_cover_the_requested_box(size):
    spacing = 0.375
    points = AutoDock4Runner._grid_points(size, spacing)
    assert points % 2 == 0
    assert points * spacing >= size
    assert points * spacing < size + 2 * spacing


@pytest.mark.parametrize("size,spacing", [(0.0, 0.375), (-1.0, 0.375), (20.0, 0.0), (float("nan"), 0.375), (20.0, float("inf"))])
def test_non_finite_or_non_positive_inputs_are_refused(size, spacing):
    with pytest.raises(ValueError):
        AutoDock4Runner._grid_points(size, spacing)
