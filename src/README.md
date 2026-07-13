# src

This directory contains the active implementation of the simulation pipeline.

## Purpose
The modules in this folder provide the core functionality for generating and shaping baseband communication signals. The current workflow includes:
- QAM modulation from binary data,
- Nyquist-based waveform synthesis,
- transmitter-side signal construction.

## Contents
- qam_modulation.py: maps bits to QAM symbols for different modulation orders.
- dac_nyquist.py: generates a Nyquist-shaped time-domain waveform from discrete symbols.
- qam_tx.py: builds a complete transmit-side signal sequence with synchronization and training components.

## Notes
This is the primary development area. New functionality should be added here, and it should remain compatible with the test suite.
