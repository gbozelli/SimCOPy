# PyCOSim RFC 0001 --- Signal, System and Component Architecture

**Status:** Draft / Spike\
**Scope:** Architecture redesign of PyCOSim\
**Audience:** PyCOSim developers\
**Source of ideas:** legacy `src/deprecated/` implementations and
`draft.py`

------------------------------------------------------------------------

## 1. Motivation

The previous PyCOSim implementation contains useful physical and
algorithmic material, but several concepts are implicit in function
arguments and array shapes.

The old code crosses multiple domains:

``` text
bits
  ↓
symbols
  ↓
electrical waveform
  ↓
optical field
  ↓
electrical waveform
  ↓
symbols
  ↓
bits
```

The new architecture should make those transitions explicit.

The main architectural goal is:

> Build a simulator in which the physical meaning of a signal, the
> capabilities of a component, and the constraints imposed by the
> selected system are explicit and validated before simulation.

This RFC is intentionally a design spike. It defines concepts and
examples before committing to a final implementation.

------------------------------------------------------------------------

# 2. Design Goals

The new architecture should:

1.  Represent digital, symbolic, electrical and optical signals.
2.  Support one or two optical polarizations.
3.  Allow multiple numerical representations of the same physical
    signal.
4.  Support both IMDD and DCS.
5.  Support WDM without making WDM-specific subclasses everywhere.
6.  Allow system-level configuration to constrain components.
7.  Validate incompatible cascades before execution.
8.  Make the signal's physical meaning explicit.
9.  Keep components composable.
10. Avoid a combinatorial explosion of subclasses.
11. Preserve the useful algorithms from the legacy implementation.
12. Make it possible to inspect a simulation configuration before
    running it.

------------------------------------------------------------------------

# 3. Non-Goals

This RFC does not yet define:

-   the final numerical implementation of every optical component;
-   exact physical equations;
-   a complete DSP framework;
-   GPU execution;
-   performance optimization;
-   a final serialization format;
-   every possible modulation format;
-   every possible optical channel model.

Those should be designed after the core abstractions are stable.

------------------------------------------------------------------------

# 4. Architectural Overview

The proposed architecture is:

``` mermaid
flowchart TD
    Simulation --> System
    Simulation --> Pipeline

    System --> Architecture
    System --> Polarization
    System --> WDM
    System --> Modulation

    Pipeline --> Transmitter
    Pipeline --> Channel
    Pipeline --> Receiver

    Transmitter --> Component
    Channel --> Component
    Receiver --> Component

    Component --> Signal
    Signal --> Domain
    Signal --> Representation
    Signal --> Polarization
```

The key relationship is:

``` text
Simulation
    │
    └── System
          │
          ├── global constraints
          │
          └── components
                  │
                  └── Signals
```

A component does not need to know everything about the simulation. It
receives the system context relevant to its validation and execution.

------------------------------------------------------------------------

# 5. Signal Is the Central Data Abstraction

The central object should be `Signal`.

A signal represents data **plus its semantic/physical meaning**.

A first conceptual API is:

``` python
class Signal:
    def __init__(
        self,
        data,
        domain,
        polarization,
        representation,
        sampling=None,
        metadata=None,
    ):
        self.data = data
        self.domain = domain
        self.polarization = polarization
        self.representation = representation
        self.sampling = sampling
        self.metadata = metadata or {}
```

Example:

``` python
signal = Signal(
    data=x,
    domain="optical",
    polarization="dual",
    representation="complex",
)
```

The important point is that the array itself does not define the
physical meaning.

------------------------------------------------------------------------

# 6. Domain

The first property of a signal is its domain.

Proposed initial domains:

``` python
from enum import Enum


class SignalDomain(Enum):
    DIGITAL = "digital"
    SYMBOL = "symbol"
    ELECTRICAL = "electrical"
    OPTICAL = "optical"
```

Examples:

``` python
bits = Signal(
    data=bits_array,
    domain=SignalDomain.DIGITAL,
    polarization=None,
    representation="binary",
)

symbols = Signal(
    data=symbol_array,
    domain=SignalDomain.SYMBOL,
    polarization=None,
    representation="complex",
)

electrical = Signal(
    data=waveform,
    domain=SignalDomain.ELECTRICAL,
    polarization="single",
    representation="real",
)

optical = Signal(
    data=field,
    domain=SignalDomain.OPTICAL,
    polarization="dual",
    representation="complex",
)
```

The domain tells us **what the signal means**, not merely what dtype it
uses.

------------------------------------------------------------------------

# 7. Domain Is Not Representation

This is one of the most important architectural rules.

A complex NumPy array is not automatically an optical signal.

Likewise, a real array with four rows is not automatically a particular
optical representation.

We therefore separate:

``` text
Domain
    What physical/semantic thing is this?

Representation
    How is that thing numerically represented?
```

Example:

``` python
signal = Signal(
    data=x,
    domain=SignalDomain.OPTICAL,
    polarization="dual",
    representation="complex",
)
```

The same physical signal may be represented as:

``` text
Complex representation:

    Ex
    Ey

where Ex and Ey are complex.
```

or:

``` text
Real representation:

    Re(Ex)
    Im(Ex)
    Re(Ey)
    Im(Ey)
```

The physical domain and polarization did not change.

------------------------------------------------------------------------

# 8. Representation

Representations should be explicit.

A first version could use:

``` python
class SignalRepresentation(Enum):
    BINARY = "binary"
    REAL = "real"
    COMPLEX = "complex"
```

Example:

``` python
optical_complex = Signal(
    data=x_complex,
    domain=SignalDomain.OPTICAL,
    polarization="dual",
    representation=SignalRepresentation.COMPLEX,
)

optical_real = Signal(
    data=x_real,
    domain=SignalDomain.OPTICAL,
    polarization="dual",
    representation=SignalRepresentation.REAL,
)
```

Both can represent the same physical optical state.

------------------------------------------------------------------------

# 9. Polarization

Polarization is a separate property.

``` python
class Polarization(Enum):
    SINGLE = "single"
    DUAL = "dual"
```

A single-polarization optical signal:

``` python
signal = Signal(
    data=ex,
    domain=SignalDomain.OPTICAL,
    polarization=Polarization.SINGLE,
    representation=SignalRepresentation.COMPLEX,
)
```

A dual-polarization signal:

``` python
signal = Signal(
    data=np.array([ex, ey]),
    domain=SignalDomain.OPTICAL,
    polarization=Polarization.DUAL,
    representation=SignalRepresentation.COMPLEX,
)
```

Conceptually:

``` text
single polarization

    Ex
    │
    └── N samples


dual polarization

    Ex
    Ey
    │
    └── N samples
```

This is preferable to making `optical_1_pol` and `optical_2_pol`
completely independent signal types.

They can still be exposed as derived properties if useful:

``` python
signal.type
# "optical_2_pol"
```

but internally the architecture should reason in terms of:

``` text
domain = optical
polarization = dual
```

------------------------------------------------------------------------

# 10. Signal Shape Must Not Define Its Meaning

Bad:

``` python
if signal.data.shape[0] == 2:
    # therefore optical dual polarization
```

Good:

``` python
if signal.domain == SignalDomain.OPTICAL:
    if signal.polarization == Polarization.DUAL:
        ...
```

The shape can be validated against the semantic metadata.

For example:

``` python
def validate_shape(signal):
    if signal.polarization == Polarization.SINGLE:
        expected_components = 1

    elif signal.polarization == Polarization.DUAL:
        expected_components = 2

    # Actual shape validation happens here.
```

This prevents silent interpretation errors.

------------------------------------------------------------------------

# 11. Signal Properties

The `Signal` class should expose useful derived information.

Example:

``` python
class Signal:
    ...

    @property
    def n_samples(self):
        return self.data.shape[-1]

    @property
    def is_complex(self):
        return self.representation == SignalRepresentation.COMPLEX

    @property
    def n_polarizations(self):
        if self.polarization == Polarization.SINGLE:
            return 1
        if self.polarization == Polarization.DUAL:
            return 2
        return 0
```

For digital and symbol signals we can additionally support:

``` python
@property
def n_symbols(self):
    ...

@property
def n_bits(self):
    ...

@property
def samples_per_symbol(self):
    ...
```

The important distinction is:

``` text
n_bits
n_symbols
n_samples
```

are different concepts.

------------------------------------------------------------------------

# 12. Sampling Metadata

The old implementation uses parameters such as samples-per-symbol and
roll-off.

These should become explicit metadata rather than undocumented
assumptions.

Example:

``` python
@dataclass
class SamplingInfo:
    symbol_rate: float | None = None
    sampling_rate: float | None = None
    samples_per_symbol: int | None = None
```

Then:

``` python
signal = Signal(
    data=x,
    domain=SignalDomain.ELECTRICAL,
    polarization=Polarization.DUAL,
    representation=SignalRepresentation.REAL,
    sampling=SamplingInfo(
        symbol_rate=32e9,
        samples_per_symbol=16,
    ),
)
```

The sampling rate can be derived:

``` python
sampling_rate = symbol_rate * samples_per_symbol
```

------------------------------------------------------------------------

# 13. Modulation Is System/Transmission Configuration

The legacy code contains QAM-specific values such as:

``` text
M
N_bits
N_sym
```

A new configuration object should make this explicit.

Example:

``` python
@dataclass
class ModulationConfig:
    name: str
    order: int

    @property
    def bits_per_symbol(self):
        return int(np.log2(self.order))
```

Then:

``` python
modulation = ModulationConfig(
    name="QAM",
    order=16,
)

print(modulation.bits_per_symbol)
# 4
```

This is preferable to scattering `M=16` throughout the pipeline.

------------------------------------------------------------------------

# 14. System

The system describes the physical architecture being simulated.

Initial proposal:

``` python
class SystemArchitecture(Enum):
    IMDD = "imdd"
    DCS = "dcs"


@dataclass
class SystemConfig:
    architecture: SystemArchitecture
    polarization: Polarization
    wdm: bool
    modulation: ModulationConfig
```

Example:

``` python
system = SystemConfig(
    architecture=SystemArchitecture.DCS,
    polarization=Polarization.DUAL,
    wdm=False,
    modulation=ModulationConfig(
        name="QAM",
        order=16,
    ),
)
```

This is a global context.

------------------------------------------------------------------------

# 15. Why System Must Be Above Components

Consider the modulator.

An IMDD modulator and a coherent modulator do not necessarily have the
same input/output semantics.

Therefore:

``` python
modulator = Modulator(system)
```

is more useful than:

``` python
modulator = Modulator(type="???")
```

The component can inspect the relevant system configuration.

Example:

``` python
class Modulator:
    def __init__(self, system):
        self.system = system

    def validate(self):
        if self.system.architecture == SystemArchitecture.IMDD:
            ...
        elif self.system.architecture == SystemArchitecture.DCS:
            ...
```

This is the beginning of configuration propagation.

------------------------------------------------------------------------

# 16. Configuration Propagation

The intended architecture is not simply "copy the System attributes into
every component".

Instead:

``` text
System
  │
  ├── defines global constraints
  │
  └── provides context
          │
          ├── Modulator
          ├── Fiber
          ├── Receiver
          └── DSP
```

A component asks the context for what it needs.

Example:

``` python
class Modulator:
    def __init__(self, context):
        self.context = context

    def validate(self):
        if self.context.system.architecture != SystemArchitecture.DCS:
            raise ConfigurationError(
                "This modulator requires a DCS system."
            )
```

This avoids duplicated configuration.

------------------------------------------------------------------------

# 17. Configuration Constraints

A system configuration should be able to forbid component
configurations.

For example:

``` text
IMDD
 ├── coherent receiver       FORBIDDEN
 ├── local oscillator        FORBIDDEN
 └── coherent optical field  FORBIDDEN

DCS
 ├── coherent receiver       ALLOWED
 ├── local oscillator        ALLOWED/REQUIRED
 └── complex optical field  REQUIRED
```

This can be formalized as validation rules.

Example:

``` python
class CoherentReceiver:
    def validate(self, system):
        if system.architecture != SystemArchitecture.DCS:
            raise ConfigurationError(
                "CoherentReceiver requires DCS."
            )
```

------------------------------------------------------------------------

# 18. Component Capabilities

Each component should declare what it accepts and produces.

A conceptual specification:

``` python
@dataclass
class SignalSpec:
    domain: SignalDomain
    polarization: Polarization | None = None
    representation: SignalRepresentation | None = None
```

Then:

``` python
class DAC:
    input_spec = SignalSpec(
        domain=SignalDomain.SYMBOL,
    )

    output_spec = SignalSpec(
        domain=SignalDomain.ELECTRICAL,
    )
```

And a coherent modulator:

``` python
class CoherentModulator:
    input_spec = SignalSpec(
        domain=SignalDomain.ELECTRICAL,
        polarization=Polarization.DUAL,
    )

    output_spec = SignalSpec(
        domain=SignalDomain.OPTICAL,
        polarization=Polarization.DUAL,
        representation=SignalRepresentation.COMPLEX,
    )
```

This creates a physical type system.

------------------------------------------------------------------------

# 19. Physical Type Checking

Suppose:

``` python
electrical_signal = Signal(
    data=x,
    domain=SignalDomain.ELECTRICAL,
    polarization=Polarization.DUAL,
    representation=SignalRepresentation.REAL,
)
```

Then:

``` python
fiber.process(electrical_signal)
```

should fail if the fiber requires optical input.

Example error:

``` text
SignalCompatibilityError:

Fiber expects:
    domain = optical

Received:
    domain = electrical
```

This is much better than allowing the numerical operation to fail later.

------------------------------------------------------------------------

# 20. Cascade Validation

Every component has:

``` text
input specification
output specification
system constraints
```

Therefore a pipeline can be checked before execution.

``` mermaid
flowchart LR
    A["Bits\nDIGITAL"] --> B["Mapper\nDIGITAL → SYMBOL"]
    B --> C["DAC\nSYMBOL → ELECTRICAL"]
    C --> D["Modulator\nELECTRICAL → OPTICAL"]
    D --> E["Fiber\nOPTICAL → OPTICAL"]
    E --> F["Receiver\nOPTICAL → ELECTRICAL"]
    F --> G["DSP\nELECTRICAL → SYMBOL"]
    G --> H["Demapper\nSYMBOL → DIGITAL"]
```

Example:

``` python
pipeline = [
    Mapper(),
    DAC(),
    CoherentModulator(),
    Fiber(),
    CoherentReceiver(),
    DSP(),
    Demapper(),
]

pipeline.validate(system)
```

The validator walks the cascade:

``` text
output(component[i])
        ↓
input(component[i+1])
```

and checks compatibility.

------------------------------------------------------------------------

# 21. Cascade Example

A simplified validator:

``` python
def validate_pipeline(pipeline, system):
    previous_output = None

    for component in pipeline:
        component.validate_system(system)

        if previous_output is not None:
            if not component.input_spec.accepts(previous_output):
                raise SignalCompatibilityError(
                    f"{component.__class__.__name__} cannot "
                    f"accept the previous signal."
                )

        previous_output = component.output_spec
```

This is only illustrative. The final implementation should probably use
a richer compatibility system.

------------------------------------------------------------------------

# 22. System Constraints + Signal Constraints

There are three levels of constraints:

``` text
Level 1
Signal
    ↓
What is this signal?

Level 2
Component
    ↓
What can this component accept/produce?

Level 3
System
    ↓
What configurations are physically valid?
```

The complete validation is:

``` text
System constraints
        ↓
Component constraints
        ↓
Signal compatibility
        ↓
Simulation execution
```

------------------------------------------------------------------------

# 23. IMDD Example

Example configuration:

``` python
system = SystemConfig(
    architecture=SystemArchitecture.IMDD,
    polarization=Polarization.SINGLE,
    wdm=False,
    modulation=ModulationConfig("OOK", 2),
)
```

A valid chain could be:

``` text
bits
  ↓
mapper
  ↓
electrical
  ↓
IMDD modulator
  ↓
optical intensity
  ↓
fiber
  ↓
photodiode
  ↓
electrical
  ↓
decision
  ↓
bits
```

Conceptual code:

``` python
pipeline = [
    BinaryMapper(),
    DAC(),
    IMDDModulator(),
    Fiber(),
    Photodiode(),
    ElectricalReceiver(),
]
```

A coherent receiver should be rejected:

``` python
pipeline = [
    BinaryMapper(),
    DAC(),
    IMDDModulator(),
    Fiber(),
    CoherentReceiver(),  # invalid
]
```

------------------------------------------------------------------------

# 24. DCS Example

Example:

``` python
system = SystemConfig(
    architecture=SystemArchitecture.DCS,
    polarization=Polarization.DUAL,
    wdm=False,
    modulation=ModulationConfig("QAM", 16),
)
```

A conceptual chain:

``` text
bits
  ↓
QAM mapper
  ↓
complex symbols
  ↓
DAC
  ↓
electrical I/Q
  ↓
coherent modulator
  ↓
complex optical field
  ↓
fiber
  ↓
coherent receiver
  ↓
electrical I/Q
  ↓
DSP
  ↓
symbols
  ↓
QAM demapper
  ↓
bits
```

Code:

``` python
pipeline = [
    QAMMapper(order=16),
    DAC(),
    CoherentModulator(),
    Fiber(),
    CoherentReceiver(),
    DSP(),
    QAMDemapper(order=16),
]
```

------------------------------------------------------------------------

# 25. Polarization Constraints

Suppose:

``` python
system = SystemConfig(
    architecture=SystemArchitecture.DCS,
    polarization=Polarization.SINGLE,
    wdm=False,
    modulation=ModulationConfig("QAM", 16),
)
```

A component requiring dual polarization should be rejected:

``` python
class DualPolCoherentModulator:
    def validate_system(self, system):
        if system.polarization != Polarization.DUAL:
            raise ConfigurationError(
                "DualPolCoherentModulator requires dual polarization."
            )
```

This is an example of a system configuration blocking another
configuration.

------------------------------------------------------------------------

# 26. Single-Pol Components Inside a Dual-Pol System

The architecture should not assume:

``` text
System = dual polarization
```

means every intermediate signal must be dual polarization.

For example:

``` text
Dual-pol system
       │
       ├── laser: carrier
       │
       ├── split into polarization branches
       │
       ├── X-pol component
       │
       ├── Y-pol component
       │
       └── recombination
```

Therefore the rule should be:

> System polarization defines the supported system configuration, while
> individual signals/components may have narrower polarization scope
> when physically meaningful.

This should be explicitly defined during implementation.

------------------------------------------------------------------------

# 27. WDM

WDM is a system-level capability:

``` python
system = SystemConfig(
    architecture=SystemArchitecture.DCS,
    polarization=Polarization.DUAL,
    wdm=True,
    modulation=ModulationConfig("QAM", 16),
)
```

A WDM signal may conceptually contain:

``` text
channel
  ├── λ1
  │    ├── X
  │    └── Y
  │
  ├── λ2
  │    ├── X
  │    └── Y
  │
  └── λ3
       ├── X
       └── Y
```

A future representation could be:

``` python
signal.data.shape
# (n_channels, n_polarizations, n_samples)
```

For example:

``` python
shape = (8, 2, 16000)
```

meaning:

``` text
8 wavelengths
2 polarizations
16000 samples/channel/polarization
```

The exact memory layout is an implementation decision, not an RFC
requirement.

------------------------------------------------------------------------

# 28. Signal Conversion

A signal may change numerical representation without changing its
physical meaning.

Example:

``` python
optical = Signal(
    data=complex_data,
    domain=SignalDomain.OPTICAL,
    polarization=Polarization.DUAL,
    representation=SignalRepresentation.COMPLEX,
)
```

Convert representation:

``` python
real_optical = optical.convert(
    representation=SignalRepresentation.REAL
)
```

The invariant is:

``` python
real_optical.domain == SignalDomain.OPTICAL
real_optical.polarization == Polarization.DUAL
```

Only the numerical representation changed.

Conceptually:

``` mermaid
flowchart LR
    A["Optical\nComplex DP"] -->|"representation conversion"| B["Optical\nReal DP"]
```

This is particularly relevant if the project later explores a
four-real-component representation.

------------------------------------------------------------------------

# 29. Hopf/Quaternion/Four-Real-Component Representations

The architecture should not commit to quaternions or Hopf coordinates
now.

Instead, it should leave room for them.

For example:

``` python
class SignalRepresentation(Enum):
    BINARY = "binary"
    REAL = "real"
    COMPLEX = "complex"
    # FUTURE:
    # HOPF = "hopf"
    # QUATERNION = "quaternion"
```

Then a future representation can be introduced without changing:

``` text
Signal.domain
Signal.polarization
System
Component
```

This is an important reason to separate domain from representation.

------------------------------------------------------------------------

# 30. Component Context

A component should receive a context rather than copying every
configuration parameter.

Example:

``` python
@dataclass
class SimulationContext:
    system: SystemConfig
```

Then:

``` python
class Fiber:
    def __init__(self, context):
        self.context = context

    def process(self, signal):
        self.validate(signal)
        ...
```

The component can access:

``` python
self.context.system.architecture
self.context.system.polarization
self.context.system.wdm
```

without duplicating those values.

------------------------------------------------------------------------

# 31. Components Should Have Local Parameters

Not everything belongs to System.

For example, fiber parameters are local:

``` python
fiber = Fiber(
    length_km=80,
    attenuation_db_per_km=0.2,
    dispersion_ps_nm_km=16.7,
)
```

The system tells the fiber:

``` text
what kind of signal/system it belongs to
```

The fiber tells itself:

``` text
what physical fiber it is
```

This distinction is important.

------------------------------------------------------------------------

# 32. Global vs Local Configuration

A useful rule:

``` text
System configuration
    = properties that define the simulation architecture

Component configuration
    = properties that define a specific physical component

Signal metadata
    = properties that define the current data being propagated
```

Example:

``` text
System:
    DCS
    dual polarization
    WDM

Modulator:
    Vpi
    bias
    extinction
    implementation model

Signal:
    optical
    dual polarization
    complex
    16 samples/symbol
```

------------------------------------------------------------------------

# 33. Example Complete Configuration

A future user-facing API could look like:

``` python
system = System(
    architecture="DCS",
    polarization="dual",
    wdm=True,
    modulation="16QAM",
)

tx = Transmitter(
    source=PRBSSource(length=2**16),
    dac=DAC(samples_per_symbol=16),
    modulator=CoherentModulator(),
)

channel = Channel(
    components=[
        Fiber(
            length_km=80,
            dispersion_ps_nm_km=16.7,
        ),
        EDFA(
            gain_db=16,
        ),
    ]
)

rx = Receiver(
    coherent=True,
    dsp=DSP(),
)

simulation = Simulation(
    system=system,
    transmitter=tx,
    channel=channel,
    receiver=rx,
)

simulation.validate()
simulation.run()
```

The important property is that:

``` python
simulation.validate()
```

can detect invalid combinations before numerical execution.

------------------------------------------------------------------------

# 34. Example of an Invalid Configuration

``` python
system = System(
    architecture="IMDD",
    polarization="single",
    wdm=False,
    modulation="OOK",
)

receiver = CoherentReceiver()

simulation = Simulation(
    system=system,
    receiver=receiver,
)

simulation.validate()
```

Expected result:

``` text
ConfigurationError:

CoherentReceiver is incompatible with
SystemArchitecture.IMDD.

Expected:
    DCS

Received:
    IMDD
```

------------------------------------------------------------------------

# 35. Compatibility Matrix

A first conceptual matrix:

  Feature                                       IMDD                     DCS
  ----------------------------- -------------------- -----------------------
  Direct intensity modulation                      ✓   possible but not core
  Complex optical field                 not required                       ✓
  Coherent receiver                                ✗                       ✓
  Local oscillator                                 ✗                       ✓
  QAM                             generally not core                       ✓
  OOK                                              ✓                possible
  Optical phase information             not required                       ✓
  Dual polarization                         possible                       ✓
  WDM                                              ✓                       ✓

This is intentionally a starting point. The final physical rules should
be reviewed before implementation.

------------------------------------------------------------------------

# 36. Component Compatibility Matrix

A first draft:

  -----------------------------------------------------------------------------------------------
  Component    Input Domain Output           Single Pol       Dual Pol           IMDD         DCS
                            Domain                                                    
  ------------ ------------ ------------ -------------- -------------- -------------- -----------
  Binary       digital      symbol                  N/A            N/A              ✓           ✓
  Mapper                                                                              

  QAM Mapper   digital      symbol                  N/A            N/A       optional           ✓

  DAC          symbol       electrical                ✓              ✓              ✓           ✓

  IMDD         electrical   optical                   ✓   configurable              ✓           ✗
  Modulator                                                                           

  Coherent     electrical   optical        configurable              ✓              ✗           ✓
  Modulator                                                                           

  Fiber        optical      optical                   ✓              ✓              ✓           ✓

  Photodiode   optical      electrical                ✓   configurable              ✓    possible

  Coherent     optical      electrical     configurable              ✓              ✗           ✓
  Receiver                                                                            

  DSP          electrical   symbol         configurable              ✓   configurable           ✓

  Demapper     symbol       digital                 N/A            N/A              ✓           ✓
  -----------------------------------------------------------------------------------------------

This table should evolve as the physical model becomes more precise.

------------------------------------------------------------------------

# 37. The Pipeline as a Typed Graph

The pipeline can be understood as a graph:

``` mermaid
flowchart LR
    B["DIGITAL"] --> M["Mapper"]
    M --> S["SYMBOL"]
    S --> DAC["DAC"]
    DAC --> E["ELECTRICAL"]
    E --> MOD["Modulator"]
    MOD --> O["OPTICAL"]
    O --> CH["Channel"]
    CH --> O2["OPTICAL"]
    O2 --> RX["Receiver"]
    RX --> E2["ELECTRICAL"]
    E2 --> DSP["DSP"]
    DSP --> S2["SYMBOL"]
    S2 --> DEM["Demapper"]
    DEM --> B2["DIGITAL"]
```

Every edge has semantic meaning.

This gives PyCOSim something similar to a type system:

``` text
component output
      ↓
must satisfy
      ↓
next component input
```

------------------------------------------------------------------------

# 38. Why This Is Better Than Raw NumPy Arrays

Without `Signal`:

``` python
fiber(signal)
```

does not tell us whether `signal` is:

``` text
bits
symbols
electrical waveform
optical field
single polarization
dual polarization
complex
real
```

With `Signal`:

``` python
fiber(signal)
```

the component can inspect:

``` python
signal.domain
signal.polarization
signal.representation
signal.sampling
```

This turns many hidden assumptions into explicit contracts.

------------------------------------------------------------------------

# 39. Relationship to the Legacy Code

The legacy material should not be thrown away conceptually.

It contains useful algorithms for:

``` text
QAM
DAC
Nyquist filtering
laser
modulation
fiber
PMD
EDFA
receiver
DSP
constellation
spectrum
```

The new architecture should wrap/reimplement those algorithms behind
explicit interfaces.

For example, instead of preserving:

``` python
DAC_Nyquist(s, SpS=16, RollOff=0.2)
```

the future interface could be:

``` python
dac = DAC(
    samples_per_symbol=16,
    roll_off=0.2,
)

output = dac.process(symbol_signal)
```

The numerical algorithm may initially be copied/adapted from the legacy
implementation.

The architecture should not depend on the old function signature.

------------------------------------------------------------------------

# 40. Proposed Simulation Lifecycle

A simulation should have explicit phases:

``` mermaid
flowchart TD
    A["Create System"] --> B["Create Components"]
    B --> C["Build Pipeline"]
    C --> D["Validate Configuration"]
    D --> E["Generate Signals"]
    E --> F["Execute Pipeline"]
    F --> G["Collect Results"]
    G --> H["Analyze"]
```

The most important phase is:

``` text
Validate Configuration
```

before numerical execution.

------------------------------------------------------------------------

# 41. Example Lifecycle

``` python
simulation = Simulation(
    system=system,
    transmitter=tx,
    channel=channel,
    receiver=rx,
)

# No numerical work yet.
simulation.validate()

# Only after successful validation:
result = simulation.run()

# Analysis:
result.plot_constellation()
result.plot_spectrum()
```

This also separates simulation from analysis.

------------------------------------------------------------------------

# 42. Analysis Should Not Be a Signal Responsibility

Rather than:

``` python
signal.constellation()
signal.spectrum()
signal.eye()
```

the first architectural preference should be analysis objects/functions:

``` python
plot_constellation(signal)
plot_spectrum(signal)
plot_eye(signal)
```

Example:

``` python
constellation = ConstellationAnalyzer(signal)
constellation.plot()
```

This keeps `Signal` focused on representing data.

------------------------------------------------------------------------

# 43. Proposed Core Package Structure

A possible future layout:

``` text
pycosim/
│
├── signal/
│   ├── signal.py
│   ├── domain.py
│   ├── representation.py
│   ├── polarization.py
│   └── sampling.py
│
├── system/
│   ├── system.py
│   ├── architecture.py
│   └── constraints.py
│
├── components/
│   ├── base.py
│   ├── source.py
│   ├── mapper.py
│   ├── dac.py
│   ├── modulator.py
│   ├── fiber.py
│   ├── edfa.py
│   ├── receiver.py
│   └── dsp.py
│
├── simulation/
│   ├── simulation.py
│   ├── pipeline.py
│   ├── validation.py
│   └── context.py
│
├── analysis/
│   ├── constellation.py
│   ├── spectrum.py
│   └── eye.py
│
└── tests/
```

This is a proposal, not a requirement.

------------------------------------------------------------------------

# 44. Suggested Base Component API

A minimal base class:

``` python
class Component:
    input_spec = None
    output_spec = None

    def validate_system(self, system):
        pass

    def validate_input(self, signal):
        if not self.input_spec.accepts(signal):
            raise SignalCompatibilityError(
                f"{type(self).__name__} cannot accept signal."
            )

    def process(self, signal):
        self.validate_input(signal)
        self.validate_system(...)
        return self._process(signal)

    def _process(self, signal):
        raise NotImplementedError
```

This creates a uniform execution model.

------------------------------------------------------------------------

# 45. Suggested SignalSpec API

``` python
@dataclass
class SignalSpec:
    domain: SignalDomain | None = None
    polarization: Polarization | None = None
    representation: SignalRepresentation | None = None

    def accepts(self, signal):
        if self.domain is not None:
            if signal.domain != self.domain:
                return False

        if self.polarization is not None:
            if signal.polarization != self.polarization:
                return False

        if self.representation is not None:
            if signal.representation != self.representation:
                return False

        return True
```

This should later become more sophisticated.

For example, a component may accept:

``` text
single OR dual polarization
```

rather than exactly one value.

------------------------------------------------------------------------

# 46. Constraints Should Be Declarative Where Possible

Instead of putting every rule inside nested `if` statements:

``` python
if system.architecture == ...:
    if system.polarization == ...:
        if system.wdm:
            ...
```

we should eventually have explicit constraints.

For example:

``` python
class CoherentReceiver(Component):

    constraints = [
        RequiresArchitecture(SystemArchitecture.DCS),
        RequiresPolarization(Polarization.DUAL),
    ]
```

Then:

``` python
component.validate_system(system)
```

evaluates those constraints.

This makes the architecture extensible.

------------------------------------------------------------------------

# 47. Configuration Inheritance vs Context

The intended meaning of "component acquires parameters from the superior
class" should be interpreted carefully.

Prefer:

``` text
System
   │
   └── Context
          │
          └── Component reads relevant configuration
```

instead of:

``` text
System
   │
   └── Component receives 25 copied attributes
```

Example:

``` python
class SimulationContext:
    system: SystemConfig

    def __init__(self, system):
        self.system = system
```

Then:

``` python
modulator = CoherentModulator(context)
```

The component has access to the context but does not own the global
configuration.

------------------------------------------------------------------------

# 48. A Concrete Future Example

The final user-facing API could eventually look approximately like:

``` python
system = System(
    architecture="DCS",
    polarization="dual",
    wdm=True,
    modulation=QAM(16),
)

simulation = Simulation(system)

simulation.transmitter(
    Source(bits=2**16),
    Mapper(),
    DAC(samples_per_symbol=16),
    CoherentModulator(),
)

simulation.channel(
    Fiber(length_km=80),
    EDFA(gain_db=16),
)

simulation.receiver(
    CoherentReceiver(),
    DSP(),
    Demapper(),
)

simulation.validate()

result = simulation.run()
```

The important design principle is that the API reads like the physical
system.

------------------------------------------------------------------------

# 49. Example of a Validation Failure

Suppose:

``` python
system = System(
    architecture="IMDD",
    polarization="single",
    wdm=False,
)
```

and the user tries:

``` python
simulation.receiver(
    CoherentReceiver()
)
```

The system should fail immediately:

``` text
ConfigurationError:
    CoherentReceiver is not valid for IMDD.

System:
    architecture = IMDD
    polarization = SINGLE

Component:
    CoherentReceiver

Required:
    architecture = DCS
```

This is preferable to discovering the problem after a numerical
simulation has already started.

------------------------------------------------------------------------

# 50. Open Questions

The following should remain open until the architecture discussion is
complete:

1.  Should `Signal` be mutable?
2.  Should `process()` always return a new `Signal`?
3.  How should multi-channel WDM data be represented?
4.  Can a component accept either single or dual polarization?
5.  Can a DP system contain single-pol intermediate signals?
6.  Should `domain` include `RF`, `analog`, or `baseband` later?
7.  Should symbols be considered a domain or a representation?
8.  Should `label` be a Signal domain or metadata?
9.  Should binary data use `bool`, integer, or packed bits?
10. Should signal conversions be automatic?
11. Should invalid conversions raise immediately?
12. How should units be represented?
13. Should physical parameters use a unit library?
14. How should component constraints be serialized?
15. Should pipeline validation be static, dynamic, or both?
16. How should stochastic components expose random seeds?
17. How should results and analysis objects be represented?

------------------------------------------------------------------------

# 51. Initial Architectural Decisions

For the first implementation spike, the following decisions are
recommended:

``` text
Decision 1
---------
Signal is the central data abstraction.

Decision 2
---------
Signal domain and numerical representation are independent.

Decision 3
---------
Polarization is a first-class Signal property.

Decision 4
---------
Single and dual polarization are first-class supported configurations.

Decision 5
---------
System defines global simulation constraints.

Decision 6
---------
Components declare input/output Signal specifications.

Decision 7
---------
Components validate against System context.

Decision 8
---------
Pipelines are validated before numerical execution.

Decision 9
---------
WDM is a System capability, with channel structure represented in Signal data/metadata.

Decision 10
---------
Avoid subclass explosion; prefer composition and configuration.

Decision 11
---------
Legacy algorithms are treated as implementation references, not as the new public architecture.

Decision 12
---------
Do not commit yet to Hopf or quaternion representations; keep the Representation abstraction extensible.
```

------------------------------------------------------------------------

# 52. Recommended First Implementation Spike

Before implementing the real optical components, build only:

``` text
Signal
SignalDomain
SignalRepresentation
Polarization
SamplingInfo
SystemConfig
SignalSpec
Component
SimulationContext
PipelineValidator
```

Then test the architecture with fake components.

Example:

``` python
source = FakeDigitalSource()

mapper = FakeMapper(
    output_domain=SignalDomain.SYMBOL
)

dac = FakeDAC()

modulator = FakeCoherentModulator()

fiber = FakeFiber()

receiver = FakeCoherentReceiver()

pipeline = [
    source,
    mapper,
    dac,
    modulator,
    fiber,
    receiver,
]

pipeline.validate(system)
```

Only after this works should the actual physical algorithms be inserted.

------------------------------------------------------------------------

# 53. Final Architecture Concept

The intended PyCOSim model can be summarized as:

``` mermaid
flowchart TD
    S["Simulation"]
    SYS["System Context"]

    SIG["Signal"]
    DOM["Domain"]
    REP["Representation"]
    POL["Polarization"]
    SAM["Sampling"]

    COMP["Component"]
    IN["Input Spec"]
    OUT["Output Spec"]
    CON["System Constraints"]

    S --> SYS
    S --> COMP

    SIG --> DOM
    SIG --> REP
    SIG --> POL
    SIG --> SAM

    COMP --> IN
    COMP --> OUT
    COMP --> CON

    SYS --> CON

    IN --> SIG
    OUT --> SIG
```

The core idea is:

> **Signals describe what is flowing. Components describe what
> transformations are possible. Systems describe which configurations
> are physically valid. Simulation coordinates the whole process.**

This gives PyCOSim a foundation where the old physical algorithms can be
rebuilt incrementally without carrying forward the old implicit
assumptions.
