#!/bin/bash

w_kessence=-0.9

# Calculate using Python
par_smg=$(python -c "print(3 * ((1 + ($w_kessence)) / 1e-07))")

echo "par_smg: $par_smg"
