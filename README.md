# sensor-cal
`sensor-cal` is a small library for calibrating magnetometers and similar devices. It performs a non-linear least squares optimization to solve for scale factor, bias, and non-orthogonality terms of your device. It works for measurements which will have a constant magnitude regardless of orientation.

## Example Usage
Collect 1000 samples from a dinky MEMS magnetometer, store them im a 1000x3 vector.

```
unsigned int max_iter = 1000;
double in_data[1000][3] = ...
```

One piece of data you'll need to know is the magnitude of your measurements, this must be roughly constant in a fixed position regardless of orientation.

```
double out_mag = 3;
```

This optimization attempts to find `A` and `b` such that `y=Ax+b` best fits your data assuming `y` has magnitude `out_mag`. You'll need to make an initial guess of `A` and `b`.

```
double A[3][3] = ...
double b[3] = ...
```

And set a max number of iterations you want the algorithm to perform.

```
double max_iter = 100;
```

Now call the `cal` function.

```
cal(in_data, in_data_len, out_mag, A, B, max_iter);
```

## Test Script
An example routine `test-cal` is provided which parses a set of data with known calibration values and runs the data through the calibration routine.
