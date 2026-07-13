# PyCOSim

PyCOSim is a Python-based simulation repository for digital communication and optical transmission experiments. The project focuses on QAM signal generation, Nyquist-shaped waveform synthesis, and validation of transmitter-side processing through automated tests.

## Repository structure
- src/: active implementation modules for the current simulation flow.
- src/deprecated/: legacy modules kept for reference and historical compatibility.
- src/tests/: unit tests for the active core modules.
- src/tests/regression/: regression tests that protect expected behavior after changes.

## Current scope
The repository currently implements:
- bit-to-QAM symbol mapping,
- waveform generation with Nyquist shaping,
- transmitter-oriented signal construction,
- automated verification through unit and regression tests.

## Development approach
The code is organized around compact, testable modules. The main development path is the active implementation under src/, while deprecated code remains available as reference material.
