from __future__ import annotations

from dataclasses import dataclass

import primer3
from Bio.Seq import Seq


class PrimerDesignError(ValueError):
    """Raised when a primer-design request cannot produce a safe result."""


@dataclass(frozen=True)
class PrimerDesignSettings:
    product_size_min: int = 100
    product_size_max: int = 800
    num_return: int = 5
    min_tm: float = 57.0
    opt_tm: float = 60.0
    max_tm: float = 63.0
    target_start: int | None = None
    target_end: int | None = None


def design_pcr_primers(
    template: str, settings: PrimerDesignSettings
) -> list[dict[str, object]]:
    """Generate forward/reverse PCR primer pairs with Primer3.

    Coordinates are zero-based half-open intervals on the provided template.
    """
    sequence = template.upper().replace(" ", "").replace("\n", "")
    invalid = sorted(set(sequence) - {"A", "C", "G", "T"})
    if invalid:
        raise PrimerDesignError(
            f"Primer design requires unambiguous A/C/G/T bases; found: {''.join(invalid)}"
        )
    if len(sequence) < settings.product_size_min:
        raise PrimerDesignError("Template is shorter than the requested minimum product size")
    if settings.product_size_min >= settings.product_size_max:
        raise PrimerDesignError("Minimum product size must be smaller than maximum product size")

    if not (50 <= settings.product_size_min < settings.product_size_max <= 10_000):
        raise PrimerDesignError("Product sizes must satisfy 50 <= minimum < maximum <= 10000")
    if not 1 <= settings.num_return <= 25:
        raise PrimerDesignError("Request between 1 and 25 primer pairs")
    if not (40 <= settings.min_tm <= settings.opt_tm <= settings.max_tm <= 80):
        raise PrimerDesignError("Temperatures must satisfy 40 <= minimum <= optimal <= maximum <= 80 °C")
    if (settings.target_start is None) != (settings.target_end is None):
        raise PrimerDesignError("Provide both target start and target end, or leave both empty")
    sequence_parameters = {"SEQUENCE_TEMPLATE": sequence}
    if settings.target_start is not None:
        start, end = settings.target_start, settings.target_end
        if not isinstance(start, int) or not isinstance(end, int) or not 0 <= start < end <= len(sequence):
            raise PrimerDesignError("Target must be a nonempty, zero-based half-open interval within the template; origin-spanning targets are not supported")
        if start < 18 or len(sequence) - end < 18:
            raise PrimerDesignError("Target needs at least 18 bases of template on each side for flanking primers")
        if end - start + 36 > settings.product_size_max:
            raise PrimerDesignError("Maximum product size cannot contain the target plus two minimum-length primers")
        sequence_parameters["SEQUENCE_TARGET"] = [start, end - start]

    result = primer3.bindings.design_primers(
        sequence_parameters,
        {
            "PRIMER_TASK": "generic",
            "PRIMER_FIRST_BASE_INDEX": 0,
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_NUM_RETURN": settings.num_return,
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_MIN_TM": settings.min_tm,
            "PRIMER_OPT_TM": settings.opt_tm,
            "PRIMER_MAX_TM": settings.max_tm,
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_PRODUCT_SIZE_RANGE": [[settings.product_size_min, settings.product_size_max]],
        },
    )
    if result.get("PRIMER_ERROR"):
        raise PrimerDesignError(str(result["PRIMER_ERROR"]))
    pairs: list[dict[str, object]] = []
    for index in range(int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))):
        left_position, left_length = result[f"PRIMER_LEFT_{index}"]
        right_five_prime, right_length = result[f"PRIMER_RIGHT_{index}"]
        pairs.append(
            {
                "pair_index": index + 1,
                "product_size": result[f"PRIMER_PAIR_{index}_PRODUCT_SIZE"],
                "left": {
                    "sequence": result[f"PRIMER_LEFT_{index}_SEQUENCE"],
                    "binding_start": left_position,
                    "binding_end": left_position + left_length,
                    "tm": result[f"PRIMER_LEFT_{index}_TM"],
                    "gc_percent": result[f"PRIMER_LEFT_{index}_GC_PERCENT"],
                    "self_any_th": result.get(f"PRIMER_LEFT_{index}_SELF_ANY_TH"),
                    "self_end_th": result.get(f"PRIMER_LEFT_{index}_SELF_END_TH"),
                },
                "right": {
                    "sequence": result[f"PRIMER_RIGHT_{index}_SEQUENCE"],
                    "binding_start": right_five_prime - right_length + 1,
                    "binding_end": right_five_prime + 1,
                    "tm": result[f"PRIMER_RIGHT_{index}_TM"],
                    "gc_percent": result[f"PRIMER_RIGHT_{index}_GC_PERCENT"],
                    "self_any_th": result.get(f"PRIMER_RIGHT_{index}_SELF_ANY_TH"),
                    "self_end_th": result.get(f"PRIMER_RIGHT_{index}_SELF_END_TH"),
                },
            }
        )
    # Independently check Primer3 coordinates against the returned oligos and scope.
    for pair in pairs:
        left, right = pair["left"], pair["right"]
        if not (0 <= left["binding_start"] < left["binding_end"] <= right["binding_start"] < right["binding_end"] <= len(sequence)):
            raise PrimerDesignError("Primer engine returned inconsistent binding coordinates")
        if left["sequence"] != sequence[left["binding_start"]:left["binding_end"]] or right["sequence"] != str(Seq(sequence[right["binding_start"]:right["binding_end"]]).reverse_complement()):
            raise PrimerDesignError("Primer engine returned oligos inconsistent with the reference")
        if pair["product_size"] != right["binding_end"] - left["binding_start"] or not settings.product_size_min <= pair["product_size"] <= settings.product_size_max:
            raise PrimerDesignError("Primer engine returned inconsistent product size")
        if settings.target_start is not None and not (left["binding_end"] <= settings.target_start < settings.target_end <= right["binding_start"]):
            raise PrimerDesignError("Primer engine returned primers that do not flank the requested target")
        pair["target"] = None if settings.target_start is None else {"start": settings.target_start, "end": settings.target_end}
        pair["design_scope"] = "linear_template_coordinates"
        pair["specificity_status"] = "not_evaluated"
    return pairs
