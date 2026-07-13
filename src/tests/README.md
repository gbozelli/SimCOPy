# tests

This directory contains automated tests for the active simulation modules.

## Purpose
The test suite validates the correctness of the core signal-processing functions and helps protect the repository from regressions during future changes.

## Contents
- test_dac_nyquist.py: checks Nyquist waveform generation behavior.
- test_qam_modulation.py: validates QAM symbol mapping.
- test_qam_tx.py: verifies the transmit-side signal generation flow.

## Notes
These tests are the primary guardrail for the current implementation. New features should be covered here whenever possible.
