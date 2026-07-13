# regression

This directory contains regression tests for long-term behavior stability.

## Purpose
Regression tests ensure that previously validated behavior remains intact after refactors, bug fixes, or feature additions.

## Contents
- test_dac_regression.py: regression coverage for DAC/Nyquist-related behavior.
- test_qam_regression.py: regression coverage for QAM mapping behavior.
- test_qam_tx_regression.py: regression coverage for the QAM transmit chain.

## Notes
These tests are intended to catch subtle changes that might not be covered by smaller unit tests.
