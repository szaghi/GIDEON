"""Derive the rows of an operator the release does not have, with the generator's own function."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'scripts'))  # the repository's scripts/
from gen_gideon import derive

# the third-order upwind-biased face value, two cells left of the face and one right; fv face: cell o spans [o-1, o]
r = derive('fv', 'face', 'point', 0, [-1, 0, 1])
print({o: str(w) for o, w in r.weights.items()}, 'exact degree', r.exact_degree, 'order', r.order)
# its mirror, one cell left and two right
r = derive('fv', 'face', 'point', 0, [0, 1, 2])
print({o: str(w) for o, w in r.weights.items()}, 'exact degree', r.exact_degree, 'order', r.order)
