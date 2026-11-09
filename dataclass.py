from dataclasses import dataclass, field
import random
import numpy as np
import math

@dataclass
class OrderedPair:
    # I'm not sure we need it, but for now it is like this
    left: int
    right: int

@dataclass
class Cohesin:
    # position is coordinates of 2 sequence positions inside the cohesin ring currently
    position: OrderedPair
    bound: bool = False

    def loop_size(self):
        # Calculates the loop size
        return self.position.right - self.position.left

    def step(self, seq, prob_block=1.0):
        # Remove the old loop end positions from the sequence
        seq.cohesin_positions.discard(self.position.left)
        seq.cohesin_positions.discard(self.position.right)
        
        # Attempt movement of right loop end one position to the right
        move_right = self.position.right + 1
        if move_right < seq.length:
            # Check for CTCF position and block probability
            blocked = (
                move_right in seq.ctcf_positions and random.random() < prob_block
            )
            if not blocked:
                self.position.right = move_right


        # Attempt movement of left loop end one posisiont to the left
        move_left = self.position.left -1
        if move_left >= 0:
            # Check for CTCF position and block probability
            blocked = (
                move_left in seq.ctcf_positions and random.random() < prob_block
            )
            if not blocked:
                self.position.left = move_left

        # Updates current positions on the sequence
        seq.cohesin_positions.add(self.position.left)
        seq.cohesin_positions.add(self.position.right)

    def bind(self, seq):
        # Attaches cohesin
        self.bound = True
        spot = seq.pick_spot()  
        self.position.left = spot
        self.position.right = spot
        seq.cohesin_positions.add(spot)

    def detach(self, seq):
        # Detaches cohesin
        self.bound = False
        seq.cohesin_positions.discard(self.position.left)
        seq.cohesin_positions.discard(self.position.right)

@dataclass 
class Sequence:
    length: int = 1000
    detachment_chance: float = 0.01
    ctcf_positions: set[int] = field(default_factory=set)
    cohesin_positions: set[int] = field(default_factory=set)

    def pick_spot(self):
        # Randomly pickes a place for attachment
        # Later can be improved by adding weights to certain locations
        spot = random.randint(0, self.length-1)
        while spot in self.cohesin_positions:
            spot = random.randint(0, self.length-1)
        return spot

@dataclass
class ContactMap:
    """
    Approximate contact map in 2D; Actual HI-C map requires 3D
    """
    sequence: Sequence
    bin_width: int
    base_prob: float = 0.1 # Probabilty of contact with no loop formation
    decay: float = 1.0 # Decay exponent due to distance between positions
    loop_multiplier: float = 2.0 # Multiplier for base probability when in the loop
    num_bins: int = field(init=False)
    total: np.ndarray = field(init=False) # Store total probabilites for all snapshots for average
    num_snapshots: int = field(init=False, default=0) # Count snapshots for average

    ### Claude helped with __post_init__
    ###     -> process of initializing values after intialization
    ###     -> choosing power law decay
    def __post_init__(self):
        # Initialize num of bins (ceil to round up for non-integers)
        self.num_bins = math.ceil(self.sequence.length/self.bin_width)
        # Initialize num_bins by num_bins matrix
        self.total = np.zeros((self.num_bins, self.num_bins))

        # Creating initial matrix (no loop)
        self._no_loop = np.zeros((self.num_bins, self.num_bins))
        for i in range(self.num_bins):
            for j in range(self.num_bins):
                dist = abs(i-j)
                # Power law decay
                self._no_loop[i][j] = self.base_prob / ((1 + dist) ** self.decay)

    def add_snapshot(self, loops):
        """
        Add contact probabilities of loops in current snapshot to total
        """
        current_snapshot_mtx = self._no_loop.copy()
        for left, right in loops:
            if right <= left:
                continue # newly bound cohesin, no loop yet
            # Find the left most and right most bin for current loop
            left_bin = left // self.bin_width
            right_bin = right //self.bin_width

            # Multiply all probabilities within the loop by multiplier
            current_snapshot_mtx[left_bin:right_bin + 1, left_bin:right_bin + 1] *= self.loop_multiplier    #fixed calculation error

        self.total += current_snapshot_mtx
        self.num_snapshots += 1

    def average(self):
        if self.num_snapshots >= 1:
            return self.total/self.num_snapshots
        else:
            return self.total

    def observed_over_expected(self):
        """
        Average map divided by the no-loop baseline. Divides out the distance baseline
        so only the effect of loops is left.
        """
        return self.average()/self._no_loop
