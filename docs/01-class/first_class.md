# First Discussion: Input Signal

Treating the input signal as a complex number vector of two units (optical and electrical signal) or a real number vector of four units (Hopf fibration), but this doesn't have an electrical or optical meaning. Quaternions can also be used, but not right now.

## Proposal 1:

```python
signal_type = [
    'electrical': 'dim -> N', 
    'optical_1_pol': 'dim -> N', 
    'optical_2_pol': 'dim -> 2N', 
    'digital': 'dim -> N * sps', 
    'label': 'dim -> N / modulation_order', 
    'binary': 'dim -> log2()'
]

pol_type = ['single', 'double']

class Signal:
    type: str  # signal_type
    pol: str   # pol_type

    def sizeof(self):
        # Calculate the size of the signal
        pass
    
    def constellation(self):
        # Each signal type has a type of constellation plot
        if self.type == 'digital':
            pass 

class System:
    WDM: bool  # True or False
    if WDM:
        Channel = 'multiple_channels'
    type: str  # 'IMDD' or 'DCS'

class Modulator:
    # This can be affected by the choice of System (IMDD or DCS).
    # What I am thinking is the use of a cascade variable.
    # With the choice of the system, we set parameters
    # and ...
    pass

class Simulation:
    # Controls the system, but maybe some configs of the system don't 
    pass
  
class Channel:
    pass
  
# The simulator can be capable of simulating IMDD and DCS
```
