# Simulation design reference

`simulation_design_reference.rds` contains only small design metadata extracted
from the user-supplied audit archive `audit_20260906_160140_13032.tar.gz`.
It contains no fitted estimates or simulation output used to populate tables.

| Source archived RData | MD5 |
| --- | --- |
| `simulate_result_pred_260717_regularC.RData` | `356f828c03d40bb1b509bb7f9be6612c` |
| `simulate_result_pred_260718_sparseC.RData` | `1ff27b0d03e42b8933209cc5d02f67ce` |

The objects were read in isolated environments. From each
`simulate_result_prediction$parameters_list[[4]]` (both multipliers equal to
one), the extraction retained the generic/sparse true C; the generic object
also supplied `Sigma_X`, `Sigma_Y`, and `VX_tilde`. Population instrument
strengths use the generic `parameters_list[c(1,3,2,4)]`, which orders the cells
as the new manuscript tables do. Source names and MD5 values are embedded in
the reference. The reference was serialized with `saveRDS(..., version = 2)`.

Script 33 uses this independent reference to check the new raw-input
calibration and truth construction. Script 32 and production data generation
do not read it and do not require either archived RData file.

Committed calibration input MD5 values:

- `paper/input/dat_1e-4.csv`: `99542ee29fef0273b637c4a11085d651`
- `paper/input/rho_mat_1e-4.csv`: `987de77f73673083a1e0d63ed30506e7`

The supplied `dat_1e-4.csv` has MD5 `37e7970c44e97a81d532b2bb4dd4bdb0`.
Only its CRLF line endings were normalized to LF in the new copy. The
correlation CSV was already LF. No numbers, column names or row order changed;
the frozen originals remain untouched.
