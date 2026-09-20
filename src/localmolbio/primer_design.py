from __future__ import annotations

from dataclasses import dataclass

import primer3


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

    result = primer3.bindings.design_primers(
        {"SEQUENCE_TEMPLATE": sequence},
        {
            "PRIMER_TASK": "generic",
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
    pairs: list[dict[str, object]] = []
    for index in range(int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))):
        left_position, left_length = result[f"PRIMER_LEFT_{index}"]
        right_three_prime, right_length = result[f"PRIMER_RIGHT_{index}"]
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
                    "binding_start": right_three_prime - right_length + 1,
                    "binding_end": right_three_prime + 1,
                    "tm": result[f"PRIMER_RIGHT_{index}_TM"],
                    "gc_percent": result[f"PRIMER_RIGHT_{index}_GC_PERCENT"],
                    "self_any_th": result.get(f"PRIMER_RIGHT_{index}_SELF_ANY_TH"),
                    "self_end_th": result.get(f"PRIMER_RIGHT_{index}_SELF_END_TH"),
                },
            }
        )
    return pairs
