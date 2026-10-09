import struct
from dataclasses import dataclass

@dataclass(frozen=True)
class BinaryFormat:
    """
    Description of a binary file format
    """

    format: str
    @property
    def size(self):
        return struct.calcsize(self.format)
    
def sink_info(feedback=True, array_length=1):
    """
    Binary format for sink_particle_info_files.
    """

    if feedback:
        main = (
            "3d 3d 3d "
            "d d d d "
            "q "
            "i i i i "
            "d d d d d d d d "
            "i "
            "200d "
            "200d "
        )
    else:
        main = (
            "3d 3d 3d "
            "d d d d "
            "q "
            "i i i i "
            "d d d d d d "
            "i "
            "200d "
            "200d "
        )

    if array_length == 1:
        extra = "d d"
    else:
        extra=f"{array_length}d {array_length}d"
    
    return BinaryFormat(main + extra)