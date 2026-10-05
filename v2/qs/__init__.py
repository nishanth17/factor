"""Exact QS/MPQS contracts, bounded collection and SIQS family primitives."""

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
from .families import FamilyStep, PolynomialFamily, family_assignments
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
from .siqs import SIQSConfig, SIQSJob

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
    "FamilyStep",
    "FilteredMatrix",
    "Polynomial",
    "PolynomialFamily",
    "PolynomialRoots",
    "PreparedRelations",
    "QSJob",
    "QSResult",
    "SIQSConfig",
    "SIQSJob",
    "SieveCollector",
    "SieveConfig",
    "SieveResult",
    "a_target",
    "build_factor_base",
    "collect_block",
    "combine_relations",
    "extract_dependency",
    "family_assignments",
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
