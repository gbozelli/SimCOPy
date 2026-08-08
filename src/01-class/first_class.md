First discussion: Input Signal

Treating the input signal as an complex number vector of two units (optical and electrical signal) or an real number vector of four units (Hopf fibration), but this dont have electrical or optical meaning. Quaternions can also be used, but not now.

Propose 1:

```python
signal_type: ['eletrical': dim->N, ['optical_1_pol': dim->N, 'optical_2_pol': dim->2N], 'digital': dim->N*sps, 'label': dim->N/modulation_order, 'binary': dim->log2()]
pol_type = [single, double]

class Signal:
  type: signal_type
  pol: pol_type

  def sizeof():
    #calculate the size of the signal
    
  def constellation():
    #each signal type have an type of constellation plot
    if type == 'digital':
      return 

class System:
  WDM: true or false
  if WDM:
    Channel = multiple_channels
  type: IMDD or DCS

class Modulator:
  # this can be affected by the choose ofSystem (IMDD or DCS)
  # what i am thinking is the use of cascade variable
  # with the choose of the system, we set parameters
  # and

class Simulation:
  # controls system, but maybe some configs of system doesnt 
  
class Channel:
  
#simulator can be capable of simulate IMDD and DCS
```
