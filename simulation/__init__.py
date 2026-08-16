"""Standardized runnable-module interface for the druggability-dossier station.

One documented command:

    python -m simulation run --input input.json --output output.json

- Reads and validates the request against schemas/input.schema.json.
- Produces the dossier via run_pipeline() (invokes the Claude Managed Agent).
- Attaches an interpretability object at dossier["interpretability"].
- Validates the dossier against schemas/output.schema.json and the interpretability
  object against schemas/interpretability.schema.json.
- Writes the dossier to --output; logs to STDERR; keeps STDOUT clean.
- Exit 0 only when both validations pass; nonzero on any failure.
"""

from .interpretability import build_interpretability
from .pipeline import (
    PipelineError,
    PipelineInvocationError,
    PipelineUnavailableError,
    run_pipeline,
)

__all__ = [
    "build_interpretability",
    "run_pipeline",
    "PipelineError",
    "PipelineUnavailableError",
    "PipelineInvocationError",
]
