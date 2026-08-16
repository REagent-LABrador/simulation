"""Standardized runnable-module interface for the druggability-dossier station.

One documented command:

    python -m simulation run --input input.json --output output.json

- Reads and validates the request against schemas/input.schema.json.
- Produces the dossier via run_pipeline(): by DEFAULT a dependency-light LOCAL
  resolver over the bundled cache in simulation/cache/ (no cloud, no Modal, no
  Paperclip, no managed agent, no API key — Python + jsonschema only). The
  managed-agent path is retained but reached only with SIMULATION_USE_AGENT=1.
- Ensures an interpretability object at dossier["interpretability"] (cached
  dossiers keep their own; others get one from build_interpretability).
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
