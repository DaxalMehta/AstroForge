from enum import IntEnum

class ParticleType(IntEnum):

    GAS = 0
    DARK_MATTER = 1
    BACKGROUND = 2
    STAR = 4
    SINK = 5

class SinkType(IntEnum):
    POPIII = 0
    POPII = 2
    BLACK_HOLE = 3