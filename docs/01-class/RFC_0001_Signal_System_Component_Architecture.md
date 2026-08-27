# PyCOSim RFC 0001 --- Signal, System and Component Architecture

**Status:** Draft / Spike\
**Version:** 2 --- merges the RFC 0002 addendum (axis semantics, composable `SignalSpec`, declarative constraints, sampling-grid compatibility, RNG discipline) directly into this document. There is no separate RFC 0002; this is now the single source of truth.\
**Scope:** Architecture redesign of PyCOSim\
**Audience:** PyCOSim developers\
**Source of ideas:** legacy `src/deprecated/` implementations, `draft.py`, and a review of the v1 draft of this RFC against `draft.py` that surfaced a concrete bug (`EDFA()` generating ASE noise from a uniform distribution instead of a Gaussian one) and several places where the v1 draft reintroduced, in its own design, the exact "shape/implicit-convention defines meaning" antipattern it set out to eliminate.

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

A review pass against `draft.py` found one concrete, currently-active
bug worth stating up front, because it motivates several decisions
below: `EDFA()` generates ASE noise using `np.random.rand` (a uniform
distribution) instead of `np.random.randn` (a Gaussian distribution).
ASE noise is physically a circular complex Gaussian process, so this
silently produces wrong noise statistics in every BER curve the script
computes. Nothing in `draft.py` could have caught this --- there are no
tests, and randomness comes from unseeded global `np.random` state.
This RFC treats "would this architecture have made that bug testable"
as a running design check, in addition to the goals in Section 2.

------------------------------------------------------------------------

## 2. Design Goals

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
13. Make stochastic components' output statistically testable
    (seeded, reproducible, per-component random streams --- see
    Section 49).

------------------------------------------------------------------------

## 3. Non-Goals

This RFC does not yet define:

-   the final numerical implementation of every optical component;
-   exact physical equations;
-   a complete DSP framework;
-   GPU execution;
-   performance optimization;
-   a final serialization format;
-   every possible modulation format;
-   every possible optical channel model;
-   a migration to an `xarray`-backed `Signal.data` (considered in
    Section 11 as a possible future direction, not committed to here).

Those should be designed after the core abstractions are stable.

------------------------------------------------------------------------

## 4. Architectural Overview

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
    Signal --> AxisSpec
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

## 5. Signal Is the Central Data Abstraction

The central object should be `Signal`.

A signal represents data **plus its semantic/physical meaning**.

``` python
class Signal:
    def __init__(
        self,
        data,
        domain,
        polarization,
        representation,
        axes: "AxisSpec | None" = None,
        sampling=None,
        metadata=None,
    ):
        self.data = data
        self.domain = domain
        self.polarization = polarization
        self.representation = representation
        self.axes = axes or AxisSpec(axes=("time",))
        self.sampling = sampling
        self.metadata = metadata or {}
        self.axes.validate_shape(data.shape)
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
physical meaning --- and, as of this version, it does not define which
array dimension means what, either. See Section 11 (`AxisSpec`) for why
this matters as soon as more than one non-time axis is involved (e.g.
polarization stacked with a WDM channel axis).

------------------------------------------------------------------------

## 6. Domain

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

## 7. Domain Is Not Representation

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

## 8. Representation

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

## 9. Polarization

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

## 10. Signal Shape Must Not Define Its Meaning

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

This rule is not limited to polarization. Section 11 generalizes it to
every axis of `signal.data`, because WDM (Section 28) turns out to
reintroduce exactly this problem on a different axis if it is not
handled the same way.

------------------------------------------------------------------------

## 11. Axis Specification

The rule in Section 10 --- "shape must not define meaning" --- applies
to every dimension of `signal.data`, not only to whether the leading
dimension is polarization. Once a second structural axis is added
(e.g. WDM channel, see Section 28), a bare `ndarray` shape such as
`(8, 2, 16000)` is again just a convention: nothing on `Signal` records
that axis 0 is "channel" and axis 1 is "polarization". A component that
receives this array has to know the convention out-of-band, which is
the same failure mode Section 10 exists to prevent.

### 11.1 Options considered

**Option A --- dict-of-signals per channel.**
Represent a multi-channel signal as `dict[int, Signal]` (channel index
→ single-channel `Signal`), each with its own `Polarization` and
`SamplingInfo`. No new abstraction is needed and channels may have
independent sampling grids, but every component now needs two code
paths (single-signal vs. dict-of-signals), pushing broadcasting logic
into every component instead of into `Signal`.

**Option B --- stacked array + explicit `AxisSpec` (recommended).**
Keep the stacked-array representation, but make axis order a validated,
first-class field on `Signal` instead of an implicit convention.

**Option C --- full `xarray`-backed `Signal.data`.**
Use `xarray.DataArray` with named dimensions and coordinates (e.g. real
wavelength values as a coordinate) instead of a bare NumPy array.
Strongest guarantees, and coordinate-aware slicing is genuinely useful
for WDM (label channels by wavelength, not index) --- at the cost of a
new runtime dependency, and it changes what "the data" is for every
component that currently expects a raw `ndarray`.

### 11.2 Recommendation

Adopt **Option B now**; keep **Option C as an explicit future
direction**, not a rejected idea (Section 3). Nothing in Option B
forecloses it, since `AxisSpec` is a thin layer that could later be
replaced by `xarray` dimension metadata without changing `Signal`'s
public surface.

``` python
from dataclasses import dataclass


@dataclass(frozen=True)
class AxisSpec:
    """Declares what each axis of Signal.data means, in order."""
    axes: tuple[str, ...]  # e.g. ("channel", "polarization", "time")

    def index_of(self, axis_name: str) -> int:
        return self.axes.index(axis_name)

    def validate_shape(self, shape: tuple[int, ...]) -> None:
        if len(shape) != len(self.axes):
            raise ShapeMismatchError(
                f"Signal declares axes {self.axes} "
                f"but data has {len(shape)} dimensions."
            )
```

A WDM signal becomes self-describing instead of convention-describing:

``` python
wdm_signal = Signal(
    data=field,  # shape (8, 2, 16000)
    domain=SignalDomain.OPTICAL,
    polarization=Polarization.DUAL,
    representation=SignalRepresentation.COMPLEX,
    axes=AxisSpec(axes=("channel", "polarization", "time")),
)

wdm_signal.axes.index_of("channel")   # 0 — no more "trust me, it's axis 0"
```

Components that are channel-agnostic (e.g. `Fiber`, applying the same
physics per channel) write `signal.axes.index_of("channel")` instead of
hardcoding `0`. A component that changes axis order (e.g. after a
transpose for FFT efficiency) expresses that explicitly by returning a
`Signal` with a different `AxisSpec`, which the pipeline validator
(Section 21) can check the same way it already checks
`domain`/`polarization`/`representation`.

------------------------------------------------------------------------

## 12. Signal Properties

The `Signal` class should expose useful derived information.

Example:

``` python
class Signal:
    ...

    @property
    def n_samples(self):
        return self.data.shape[self.axes.index_of("time")]

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

## 13. Sampling Metadata

The old implementation uses parameters such as samples-per-symbol and
roll-off.

These should become explicit metadata rather than undocumented
assumptions, and it should be possible to check whether two signals
share a grid rather than assuming they do --- see Section 48 for why
this matters as soon as a component combines two signals.

``` python
@dataclass
class SamplingInfo:
    symbol_rate: float | None = None
    sampling_rate: float | None = None
    samples_per_symbol: int | None = None

    def compatible_with(self, other: "SamplingInfo", rtol: float = 1e-9) -> bool:
        if self.sampling_rate is None or other.sampling_rate is None:
            return False
        return abs(self.sampling_rate - other.sampling_rate) <= rtol * self.sampling_rate
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

## 14. Modulation Is System/Transmission Configuration

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

## 15. System

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

## 16. Why System Must Be Above Components

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

## 17. Configuration Propagation

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

## 18. Configuration Constraints

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

This can be formalized as validation rules. See Section 45--47 for the
single declarative mechanism used to express these rules.

------------------------------------------------------------------------

## 19. Component Capabilities

Each component should declare what it accepts and produces.

A conceptual specification (the full, composable version used from
implementation onward is given in Section 46):

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

## 20. Physical Type Checking

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

## 21. Cascade Validation

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

## 22. Cascade Example

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

This is only illustrative. The final implementation uses the
`SignalSpec.accepts()` given in Section 46, which supports "accepts
single OR dual polarization" rather than exact match only.

------------------------------------------------------------------------

## 23. System Constraints + Signal Constraints

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

## 24. IMDD Example

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

## 25. DCS Example

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

## 26. Polarization Constraints

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
configuration. As of Section 45, this rule is expressed declaratively
(`RequiresPolarization`) rather than as a hand-written override.

------------------------------------------------------------------------

## 27. Single-Pol Components Inside a Dual-Pol System

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

## 28. WDM

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

The stacked-array representation is:

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

Unlike the earlier draft of this RFC, this shape is not treated as a
convention that components must know out-of-band: it is declared via
`AxisSpec` (Section 11) and carried on the `Signal` instance itself:

``` python
wdm_signal = Signal(
    data=field,
    domain=SignalDomain.OPTICAL,
    polarization=Polarization.DUAL,
    representation=SignalRepresentation.COMPLEX,
    axes=AxisSpec(axes=("channel", "polarization", "time")),
)
```

The exact memory layout is still an implementation decision, but which
axis means what is no longer one.

------------------------------------------------------------------------

## 29. Signal Conversion

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

## 30. Hopf/Quaternion/Four-Real-Component Representations

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

## 31. Component Context

A component should receive a context rather than copying every
configuration parameter.

``` python
@dataclass
class SimulationContext:
    system: SystemConfig
    seed: int

    def __post_init__(self):
        self._seed_sequence = np.random.SeedSequence(self.seed)
        self._issued: dict[str, "np.random.Generator"] = {}

    def rng_for(self, component_name: str) -> "np.random.Generator":
        """Each component gets its own independent, reproducible stream."""
        if component_name not in self._issued:
            child_seed = self._seed_sequence.spawn(1)[0]
            self._issued[component_name] = np.random.default_rng(child_seed)
        return self._issued[component_name]
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

without duplicating those values, and, for stochastic components, a
private, reproducible random stream via `self.context.rng_for(...)`
instead of module-level `np.random` state (see Section 49).

------------------------------------------------------------------------

## 32. Components Should Have Local Parameters

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

## 33. Global vs Local Configuration

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

`Signal.metadata` (Section 5) is a free-form escape hatch and should
stay one deliberately: anything that needs to be validated or relied
upon by components --- axis layout, sampling grid, RNG seed --- belongs
in a typed field (`AxisSpec`, `SamplingInfo`, `SimulationContext`), not
in `metadata`. Otherwise `metadata` quietly becomes a second
implicit-shape problem of the kind Section 10 exists to prevent.

------------------------------------------------------------------------

## 34. Example Complete Configuration

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

## 35. Example of an Invalid Configuration

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

## 36. Compatibility Matrix

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

## 37. Component Compatibility Matrix

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

Note that several components (`DAC`, `Fiber`, `Photodiode`) accept
**either** single **or** dual polarization. Section 46's `SignalSpec`
is written to express exactly this ("accepts A or B"), which the
Section 19 sketch could not.

This table should evolve as the physical model becomes more precise.

------------------------------------------------------------------------

## 38. The Pipeline as a Typed Graph

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

## 39. Why This Is Better Than Raw NumPy Arrays

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
signal.axes
```

This turns many hidden assumptions into explicit contracts.

------------------------------------------------------------------------

## 40. Relationship to the Legacy Code

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
implementation, but not verbatim: the review pass found `fiber()`
duplicated identically twice, `QAM_receiver_DP` and `DPQAM_receiver`
each defined twice with the later definition silently shadowing the
earlier one, and the `EDFA()` ASE-noise bug described in Section 1.
These should be fixed --- not carried forward --- as each algorithm is
ported behind the new interfaces, and the Gray-coding lookup tables
(currently ~230 lines of hand-written `if/elif` branches per modulation
order) should become a generator function rather than literal tables.

The architecture should not depend on the old function signature.

------------------------------------------------------------------------

## 41. Proposed Simulation Lifecycle

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

## 42. Example Lifecycle

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

## 43. Analysis Should Not Be a Signal Responsibility

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

This keeps `Signal` focused on representing data. It also keeps
component functions free of the side effects found throughout
`draft.py` (`matplotlib` calls embedded in physics functions via
`plot_flag` arguments), which makes headless/batch execution and
testing harder than necessary.

------------------------------------------------------------------------

## 44. Proposed Core Package Structure

A possible future layout:

``` text
pycosim/
│
├── signal/
│   ├── signal.py
│   ├── domain.py
│   ├── representation.py
│   ├── polarization.py
│   ├── axes.py
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

## 45. Suggested Base Component API

There is exactly one validation mechanism, used by every component:
a declarative `constraints` list, evaluated once by the base class.
`Component.validate_system` is never overridden by a subclass ---
subclasses only declare data (`constraints`, `input_spec`,
`output_spec`). This is what makes a simulation configuration
inspectable (Design Goal 12): the constraint list can be walked and
rendered into an error report (Section 52) without parsing subclass
source code, which an imperative `if`-based override cannot support.

``` python
class Constraint:
    def check(self, system: "SystemConfig") -> None:
        """Raise ConfigurationError if not satisfied."""
        raise NotImplementedError


@dataclass
class RequiresArchitecture(Constraint):
    architecture: "SystemArchitecture"

    def check(self, system):
        if system.architecture != self.architecture:
            raise ConfigurationError(
                f"{self.architecture} required, got {system.architecture}."
            )


@dataclass
class RequiresPolarization(Constraint):
    polarization: "Polarization"

    def check(self, system):
        if system.polarization != self.polarization:
            raise ConfigurationError(
                f"{self.polarization} required, got {system.polarization}."
            )
```

``` python
class Component:
    input_spec: "SignalSpec" = SignalSpec()
    output_spec: "SignalSpec" = SignalSpec()
    constraints: list[Constraint] = []

    def validate_system(self, system: "SystemConfig") -> None:
        for constraint in self.constraints:
            constraint.check(system)

    def validate_input(self, signal):
        if not self.input_spec.accepts(signal):
            raise SignalCompatibilityError(
                f"{type(self).__name__} cannot accept signal."
            )

    def process(self, signal, system):
        self.validate_input(signal)
        self.validate_system(system)
        return self._process(signal)

    def _process(self, signal):
        raise NotImplementedError
```

A concrete component declares only data:

``` python
class CoherentReceiver(Component):
    constraints = [
        RequiresArchitecture(SystemArchitecture.DCS),
        RequiresPolarization(Polarization.DUAL),
    ]
```

This creates a uniform execution model, and gives the error report in
Section 52 for free --- walk `component.constraints`, render each one
--- instead of a bespoke message per subclass.

------------------------------------------------------------------------

## 46. Suggested SignalSpec API

`SignalSpec` fields accept `None` (any value), an exact value, or a
`frozenset` of allowed values --- never exact-match-only. This is
required from the first implementation, not deferred: the Section 37
compatibility table already declares `Fiber`, `DAC`, and `Photodiode`
as accepting single **or** dual polarization, so an exact-match-only
`accepts()` cannot express the components that appear in every pipeline
example in this document.

``` python
from typing import Union

PolarizationConstraint = Union["Polarization", "frozenset[Polarization]", None]


@dataclass
class SignalSpec:
    domain: "SignalDomain | frozenset[SignalDomain] | None" = None
    polarization: PolarizationConstraint = None
    representation: "SignalRepresentation | frozenset[SignalRepresentation] | None" = None

    def accepts(self, signal) -> bool:
        return (
            self._field_ok(self.domain, signal.domain)
            and self._field_ok(self.polarization, signal.polarization)
            and self._field_ok(self.representation, signal.representation)
        )

    @staticmethod
    def _field_ok(constraint, value) -> bool:
        if constraint is None:
            return True
        if isinstance(constraint, frozenset):
            return value in constraint
        return value == constraint
```

`Fiber` now matches the compatibility table it is supposed to satisfy:

``` python
class Fiber(Component):
    input_spec = SignalSpec(
        domain=SignalDomain.OPTICAL,
        polarization=frozenset({Polarization.SINGLE, Polarization.DUAL}),
    )
    output_spec = SignalSpec(
        domain=SignalDomain.OPTICAL,
        polarization=frozenset({Polarization.SINGLE, Polarization.DUAL}),
    )
```

------------------------------------------------------------------------

## 47. Constraints Should Be Declarative

Instead of nested `if` statements scattered per subclass:

``` python
if system.architecture == ...:
    if system.polarization == ...:
        if system.wdm:
            ...
```

every component declares its rules as data, evaluated generically by
the base class described in Section 45:

``` python
class CoherentReceiver(Component):
    constraints = [
        RequiresArchitecture(SystemArchitecture.DCS),
        RequiresPolarization(Polarization.DUAL),
    ]
```

This is the only validation path (Section 45) --- there is no separate
imperative override to keep in sync with it --- which makes the
architecture extensible: a new rule is a new `Constraint` subclass, not
a new `if` branch buried in a component that also does physics.

------------------------------------------------------------------------

## 48. Sampling Grid Compatibility

`SamplingInfo` (Section 13) can report whether two signals share a
grid, but nothing yet requires a component to check this before
combining signals. This matters for any component that mixes two
inputs --- a coherent hybrid combining a signal with a local
oscillator, a WDM multiplexer, a polarization combiner --- since
nothing in `draft.py` prevents combining signals sampled on different
grids; it "works" only because every example in the notebook happens
to reuse the same `ts`/`SpS` throughout.

``` python
@dataclass
class RequiresSameGrid(Constraint):
    """Component-level constraint: all input signals must share a time grid."""

    def check_signals(self, signals: list["Signal"]) -> None:
        reference = signals[0].sampling
        for other in signals[1:]:
            if not reference.compatible_with(other.sampling):
                raise SignalCompatibilityError(
                    "Inputs do not share a common sampling grid: "
                    f"{reference} vs {other.sampling}."
                )
```

Any component combining multiple signals declares `RequiresSameGrid()`
alongside its other `constraints` (Section 45) and calls
`check_signals()` on its inputs before combining them.

------------------------------------------------------------------------

## 49. Random Number Generation Discipline

This is treated as a Day-1 decision rather than deferred (compare the
open question list in Section 53), because it is the decision that
would have made the Section 1 bug (`EDFA()` drawing ASE noise from a
uniform instead of a Gaussian distribution) catchable in a unit test
instead of invisible. Global, unseeded `np.random` state --- used
throughout `draft.py` --- cannot be asserted against; a per-component,
seeded, injectable generator can.

`SimulationContext` (Section 31) owns RNG derivation. Components never
read or write module-level `np.random` state:

``` python
class EDFA(Component):
    def __init__(self, context, gain_db, noise_figure_db):
        self.rng = context.rng_for(f"EDFA:{id(self)}")
        self.gain_db = gain_db
        self.noise_figure_db = noise_figure_db

    def _process(self, signal):
        ...
        noise = (
            self.rng.standard_normal(n_samples)
            + 1j * self.rng.standard_normal(n_samples)
        ) * noise_std
        ...
```

With this in place, a unit test can seed the context, run `EDFA` on a
zero-input signal, and assert the output noise samples pass a normality
check (e.g. `scipy.stats.normaltest`) --- something that cannot be
written against `np.random.rand`-based noise, or reliably against
unseeded global state, regardless of which distribution is used.

Two consequences follow and are decided here, not deferred:

- Seeds are always passed explicitly to `SimulationContext`; there is
  no implicit "random by default" mode in library code (a top-level
  CLI/script convenience wrapper may choose one, but that is outside
  this RFC's scope).
- `RequiresSeededContext` is added as a standard `Constraint`
  (Section 45) so stochastic components can declare that dependency
  the same way they declare architecture or polarization requirements.

------------------------------------------------------------------------

## 50. Configuration Inheritance vs Context

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

`SimulationContext` (Section 31) is this context object, and now also
owns RNG derivation (Section 49) for the same reason: components read
what they need from a single shared object rather than each holding
its own copy of simulation-wide state.

------------------------------------------------------------------------

## 51. A Concrete Future Example

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

## 52. Example of a Validation Failure

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

## 53. Open Questions

The following should remain open until the architecture discussion is
complete:

1.  Should `Signal` be mutable?
2.  Should `process()` always return a new `Signal`?
3.  Should `AxisSpec` allow *named* channel coordinates (e.g. actual
    wavelength values), or stay purely positional-with-labels until/
    unless the `xarray` migration (Section 11) happens?
4.  Can a DP system contain single-pol intermediate signals?
5.  Should `domain` include `RF`, `analog`, or `baseband` later?
6.  Should symbols be considered a domain or a representation?
7.  Should `label` be a Signal domain or metadata?
8.  Should binary data use `bool`, integer, or packed bits?
9.  Should signal conversions be automatic?
10. Should invalid conversions raise immediately?
11. How should units be represented?
12. Should physical parameters use a unit library?
13. How should component constraints be serialized?
14. Should `rng_for()` key by component *instance* (current proposal,
    Section 49) or by component *class + position in pipeline*, to
    keep streams stable across pipeline edits that don't touch the
    component itself?
15. Should `RequiresSameGrid` (Section 48) allow a per-component
    tolerance, or should the global `rtol` in
    `SamplingInfo.compatible_with` be the single source of truth?
16. How should results and analysis objects be represented?

Resolved by this version (previously open, see the original v1 list):

-   *"Can a component accept either single or dual polarization?"* ---
    resolved, Section 46 (`SignalSpec` supports `frozenset`).
-   *"How should stochastic components expose random seeds?"* ---
    resolved, Section 49 (`SimulationContext.rng_for`).
-   *"How should multi-channel WDM data be represented?"* --- resolved
    for axis semantics, Section 11 (`AxisSpec`); memory layout details
    remain open.
-   *"Should pipeline validation be static, dynamic, or both?"* ---
    effectively static-before-run via the Section 45 constraint list;
    dynamic per-signal checks still happen in `validate_input`/
    `process`, so both remain in play.

------------------------------------------------------------------------

## 54. Initial Architectural Decisions

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
WDM is a System capability. Channel structure is represented via an
explicit AxisSpec on Signal, not a positional shape convention.

Decision 10
---------
Avoid subclass explosion; prefer composition and configuration.

Decision 11
---------
Legacy algorithms are treated as implementation references, not as the
new public architecture, and known legacy bugs (duplicate function
definitions, the EDFA noise-distribution bug) are fixed during porting,
not carried forward.

Decision 12
---------
Do not commit yet to Hopf or quaternion representations; keep the
Representation abstraction extensible.

Decision 13
---------
Signal carries an explicit AxisSpec; axis order is declared, not
assumed.

Decision 14
---------
SignalSpec fields accept None (any), an exact value, or a frozenset of
allowed values — never exact-match-only.

Decision 15
---------
Component.validate_system is implemented once, in the base class, and
iterates a declarative constraints list. No subclass overrides it.

Decision 16
---------
Components that combine multiple signals declare RequiresSameGrid and
check sampling-grid compatibility before combining.

Decision 17
---------
SimulationContext owns RNG derivation (numpy SeedSequence.spawn per
component). No component reads or writes global np.random state.

Decision 18
---------
An xarray-backed Signal.data remains a possible future direction and is
explicitly not decided here (see Section 3, Section 11).
```

------------------------------------------------------------------------

## 55. Recommended First Implementation Spike

Before implementing the real optical components, build only:

``` text
AxisSpec
Signal (axis-aware)
SignalDomain
SignalRepresentation
Polarization
SamplingInfo (+ compatible_with)
SignalSpec (composable — supports sets, not just exact match)
Constraint, RequiresArchitecture, RequiresPolarization, RequiresSameGrid
SystemConfig
Component (base class owns validate_system via constraints)
SimulationContext (+ rng_for)
PipelineValidator
```

Then test the architecture with fake components:

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

Because `SignalSpec` now supports sets (Section 46), the fake
`Fiber`/`DAC` can legitimately declare
`polarization=frozenset({SINGLE, DUAL})`; because `Component.
validate_system` lives only in the base class (Section 45), no fake
needs its own override; and a fake stochastic component can exercise
`context.rng_for(...)` (Section 49) to prove seeding works end-to-end
before any real physics is ported.

Only after this works should the actual physical algorithms --- ported
per Section 40, with the known legacy bugs fixed rather than
reproduced --- be inserted behind these interfaces.

------------------------------------------------------------------------

## 56. Final Architecture Concept

The intended PyCOSim model can be summarized as:

``` mermaid
flowchart TD
    S["Simulation"]
    SYS["System Context"]
    RNG["RNG per Component"]

    SIG["Signal"]
    DOM["Domain"]
    REP["Representation"]
    POL["Polarization"]
    SAM["Sampling"]
    AX["AxisSpec"]

    COMP["Component"]
    IN["Input Spec"]
    OUT["Output Spec"]
    CON["Constraints"]

    S --> SYS
    S --> COMP
    SYS --> RNG

    SIG --> DOM
    SIG --> REP
    SIG --> POL
    SIG --> SAM
    SIG --> AX

    COMP --> IN
    COMP --> OUT
    COMP --> CON

    SYS --> CON

    IN --> SIG
    OUT --> SIG
```

The core idea is:

> **Signals describe what is flowing, including which axis means what.
> Components describe what transformations are possible, and declare
> their rules as data rather than code. Systems describe which
> configurations are physically valid. Simulation coordinates the whole
> process, including reproducible randomness.**

This gives PyCOSim a foundation where the old physical algorithms can be
rebuilt incrementally without carrying forward the old implicit
assumptions --- or its one confirmed bug.
