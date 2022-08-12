"""HiFi Assembly

Description of HiFi Assembly

Usage Example
    ...example
"""
# Standard Library


# FUNCTIONS -------------------------------------------------------------------


# CLASSES ---------------------------------------------------------------------
class Fragment:
    """placeholder"""

    def __init__(self):
        fragment_id=0
        project_id = 0
        fragment_name = 0
        method = 0
        source = None
        primer_1 = None
        primer_2 = None
        bp = None
        notes = None
        sequence = None


class Construct:
    """placeholder"""

    def __init__(self, *args):
        construct_id=None
        project_id=None
        construct_name=None
        method=None
        bp=None
        map_created=None
        designed_by=None
        location=None
        sequence_verified=None
        cloned_by=None
        date_banked=None
        banked_by=None
        strain=None
        notes=None
        parts=args


class HiFiAssembly:
    """placeholder"""

    def __init__(self):
        reaction_volume=0
        vector_concentration=0
        construct=None

    def parts(self):
        for part in self.construct.parts:
            pass
        return
