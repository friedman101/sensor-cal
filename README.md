# sensor-cal
`sensor-cal` is a small library for calibrating magnetometers and similar devices. It performs a non-linear least squares optimization to solve for scale factor, bias, and non-orthogonality terms of your device. It works for measurements which will have a constant magnitude regardless of orientation.

The assumption is the relationship between the engineering unit output of your device (say, magnetic field in Guass) and actual output (say, volts) can be described by a linear transform of the type

```
y = Ax + b
```

Where `A` is a matrix that contains scale factor and non-orthogonality terms

```
A = [SFx 0 0; yx SFy 0; zx zy SFz]
```

and `b` are bias values
```
b = [bx; by; bz]
```

Why does `A` contain 3 cells of zeros? As this is an entirely referenceless calibration the axes can only be misaligned relative to one another. The x-axis is the reference axis and thus cannot be misaligned by definition. The y-axis is only partially defined as it is expected to by 90° from x. Any deviation from this results in a non-zero `yz`. The z-axis is then fully defined as `cross(x,y)` and any deviation from this results in a non-zero `zx` and/or non-zero `zy`.

## Building
The only dependancies are libgsl and cmake. Use brew, apt-get, yum, etc to get those. Then just clone and build

```
git clone git@github.com:friedman101/sensor-cal.git
cd sensor-cal
cmake .
make
```

## Example Usage
Collect 1000 samples from a dinky MEMS magnetometer, store them in a 1000x3 vector.

```
unsigned int x_len = 1000;
double x[1000][3] = ...
```

One piece of data you'll need to know is the magnitude of your measurements, this must be roughly constant in a fixed position regardless of orientation.

```
double y_mag = 3;
```

This optimization attempts to find `A` and `b` such that `y=Ax+b` best fits your data assuming `y` has magnitude `y_mag`. You'll need to make an initial guess of `A` and `b`.

```
double A[3][3] = ...
double b[3] = ...
```

Set a max number of iterations you want the algorithm to perform.

```
double max_iter = 100;
```

And a tolerance for when to stop iterating. Once the algorithm's solution is improving by less than `tol*100` percent it will stop iterating.

```
double tol = 1e-4;
```


Now call the `cal` function.

```
cal(x, x_len, y_mag, A, b, max_iter, tol);
```

## Test Utility
An example routine `test-cal` is provided which parses a set of data with known calibration values and runs the data through the calibration routine.

To generate the test vectors you'll need to run `genTestData` using Octave.

Running `test-cal` on my machine produces the following output

```
before cal:
A = [1.500000 0.000000 0.000000;
 2.500000 3.500000 0.000000;
 4.500000 5.500000 6.500000]
b = [7.500000; 8.500000; 9.500000]
took 273.86 ms to optimize 100000 samples
done in 9 iterations
after cal:
A = [0.999972 0.000000 0.000000;
 1.999765 2.999738 0.000000;
 4.000188 5.000156 5.999929]
b = [6.999744; 7.998878; 9.001096]
```
