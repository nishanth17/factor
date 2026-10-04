"""Exact QS/MPQS contracts and bounded collection; SIQS comes in P3.4."""

from .extraction import (
    CongruenceResult,
    DependencyExtractor,
    PreparedRelations,
    extract_dependency,
    prepare_relations,
)
from .factor_base import (
    FactorBase,
    FactorBaseBuild,
    FactorBaseEntry,
    build_factor_base,
    modular_square_roots,
)
from .linear_algebra import (
    DependencySolver,
    FilteredMatrix,
    filter_matrix,
    verify_dependency,
)
from .pipeline import QSJob, QSResult
from .polynomial import (
    Polynomial,
    PolynomialRoots,
    a_target,
    mpqs_polynomial,
    polynomial_roots,
    qs_polynomial,
)
from .reference_collector import CollectionResult, collect_block
from .relations import (
    AtomicRelation,
    CombinationResult,
    CombinedRelation,
    combine_relations,
    parity_bits,
    verify_atomic,
    verify_combined,
)
from .sieve_collector import SieveCollector, SieveConfig, SieveResult

__all__ = [
    "AtomicRelation",
    "CollectionResult",
    "CombinationResult",
    "CombinedRelation",
    "CongruenceResult",
    "DependencyExtractor",
    "DependencySolver",
    "FactorBase",
    "FactorBaseBuild",
    "FactorBaseEntry",
    "FilteredMatrix",
    "Polynomial",
    "PolynomialRoots",
    "PreparedRelations",
    "QSJob",
    "QSResult",
    "SieveCollector",
    "SieveConfig",
    "SieveResult",
    "a_target",
    "build_factor_base",
    "collect_block",
    "combine_relations",
    "extract_dependency",
    "filter_matrix",
    "modular_square_roots",
    "mpqs_polynomial",
    "parity_bits",
    "polynomial_roots",
    "prepare_relations",
    "qs_polynomial",
    "verify_atomic",
    "verify_combined",
    "verify_dependency",
]
