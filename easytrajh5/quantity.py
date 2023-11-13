import numpy as np
from parmed import unit
from parmed.unit import *


def compile_unit(unit_repr):
    """
    This is a hack to allow the generation of a unit from a string representation.
     It uses eval() to reconstruct the unit, and so requires all available units
     in the immediate scope. This has been pushed into this module to avoid the
     clobbering of the namespace from the unit module, and to isolate the
     usage of eval
    """
    return eval(unit_repr)


def get_dict_from_quantity(o):
    value = o._value
    if isinstance(value, np.ndarray):
        value = value.tolist()
    elif isinstance(value, list):
        value = np.array(value).tolist()
    result = {"type": "quantity", "value": value, "unit": str(o.unit)}
    if not hasattr(unit, result["unit"]):
        result["unit_repr"] = repr(o.unit)
    return result


def get_quantity_from_dict(o):
    if hasattr(unit, o["unit"]):
        return o["value"] * getattr(unit, o["unit"])
    elif o.get("unit_repr"):
        return o["value"] * compile_unit(o["unit_repr"])
    return o
