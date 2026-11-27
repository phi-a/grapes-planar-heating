from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Vehicle:
    name: str
    mass: float  # kg
    ref_area: float  # m^2
    cd: float
    nose_radius: float  # m
    lift_to_drag: float = 0.0

    @property
    def ballistic_coefficient(self) -> float:
        """m / (Cd A), kg/m^2."""
        return self.mass / (self.cd * self.ref_area)
