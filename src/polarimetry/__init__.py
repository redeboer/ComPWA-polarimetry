from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import sympy as sp
from ampform.sympy import PoolSum
from sympy.physics.matrices import msigma

from polarimetry.spin import create_spin_range

if TYPE_CHECKING:
    from .amplitude import DalitzPlotDecompositionBuilder


def formulate_polarimetry(
    builder: DalitzPlotDecompositionBuilder, reference_subsystem: Literal[1, 2, 3] = 1
) -> tuple[PoolSum, PoolSum, PoolSum]:
    half = sp.Rational(1, 2)
    if builder.decay.initial_state.spin != half:
        msg = (
            "Can only formulate polarimetry for an initial state with spin 1/2, but"
            f" got {builder.decay.initial_state.spin}"
        )
        raise ValueError(msg)
    model = builder.formulate(reference_subsystem)
    λ0, λ0_prime = sp.symbols(R"lambda \lambda^{\prime}", rational=True)
    λ = {
        sp.Symbol(f"lambda{i}", rational=True): create_spin_range(state.spin)
        for i, state in builder.decay.final_state.items()
    }
    ref = reference_subsystem
    return tuple(
        PoolSum(
            builder.formulate_aligned_amplitude(λ0, *λ, ref)[0].conjugate()
            * pauli_matrix[_to_index(λ0), _to_index(λ0_prime)]
            * builder.formulate_aligned_amplitude(λ0_prime, *λ, ref)[0],
            (λ0, [-half, +half]),
            (λ0_prime, [-half, +half]),
            *λ.items(),
        ).cleanup()
        / model.intensity
        for pauli_matrix in map(msigma, [1, 2, 3])
    )


def _to_index(helicity):
    """Symbolic conversion of half-value helicities to Pauli matrix indices."""
    return sp.Piecewise(
        (1, sp.LessThan(helicity, 0)),
        (0, True),
    )
