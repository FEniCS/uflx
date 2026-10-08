# Benchmarks

A first run of the benchmarks of both pipelines, measured on 2026-10-07 at commit
`9e8efaf`.

## Machine

- CPU: two AMD EPYC 7742 sockets (Zen 2), with 64 cores each and one thread per
  core. The base clock is 2.25 GHz and the maximum boost is 3.4 GHz. Each core does
  two 256-bit FMAs per cycle, which is 16 double precision flops.
- Memory: 2 TB of DDR4 LRDIMM over 16 channels, in two NUMA nodes.
- GPU: NVIDIA Tesla V100-PCIE-32GB, of which the benchmarks use one. It has 80 SMs
  at up to 1380 MHz and 900 GB/s of HBM2.
- Software: Debian 13, LLVM and MLIR 22.1.2 with the patches in the README, CUDA
  12.4 and driver 550.163.01.

The peaks used below are as follows.

| | FP64 peak | Memory bandwidth |
|---|---|---|
| V100 | 7.07 TFLOP/s | 798 GB/s, measured with a copy kernel |
| 128 CPU cores | 4.61 TFLOP/s at the base clock | 284 GB/s, measured with a triad |
| 64 CPU cores | 2.30 TFLOP/s at the base clock | 284 GB/s, measured with a triad |

The triad bandwidth counts the write-allocate traffic of the stored array. Counting
only the three arrays it is 213 GB/s. Under load the cores run above their base
clock, so the CPU's real compute peak is higher than in the table, by up to 50%.

## Setup

- Forms: mass, Laplace and Helmholtz, on hexahedra, with GLL-warped Lagrange
  elements of degree 1 to 7.
- Quadrature of degree 2k + 1, which has (k + 1)³ points, as many as there are dofs.
- Straight-sided cells, with geometry of degree 1.
- `pre` reads geometry data precomputed at each quadrature point. `inline` reads the
  cell's vertices and computes the geometry in the kernel.
- The action acts on E-vectors, so there is no gather or scatter, except in the
  section on gather and scatter.
- GPU: 3.2 × 10⁷ dofs, with the best of 1, 2, 4, 8 and 16 cells per block. The
  largest block size that compiles at each degree fails to launch, because its
  static shared memory is over 48 KB, so it is left out.
- CPU: 6.4 × 10⁷ dofs, with the best of 4 and 8 cells per vector block. OpenMP
  threads are pinned and spread over both sockets, and each thread first touches
  its own data.
- Each time is the mean over at least one second of calls. Every kernel passed its
  check against the NumPy reference.

## Results

Throughput in GDoF/s on the GPU.

| Q | Mass, pre | Laplace, pre | Helmholtz, pre | Mass, inline | Laplace, inline | Helmholtz, inline |
|---|---|---|---|---|---|---|
| 1 | 33.14 | 11.60 | 10.34 | 18.11 | 14.20 | 11.98 |
| 2 | 31.27 | 10.82 | 9.60 | 23.36 | 11.51 | 10.26 |
| 3 | 31.54 | 11.18 | 9.05 | 26.67 | 12.74 | 9.61 |
| 4 | 32.64 | 11.19 | 9.36 | 28.42 | 12.41 | 10.92 |
| 5 | 32.89 | 11.12 | 8.78 | 21.25 | 7.67 | 6.86 |
| 6 | 32.18 | 10.40 | 8.73 | 26.57 | 9.54 | 8.34 |
| 7 | 32.73 | 10.87 | 9.09 | 30.95 | 11.55 | 9.32 |

Throughput in GDoF/s on 128 CPU cores.

| Q | Mass, pre | Laplace, pre | Helmholtz, pre | Mass, inline | Laplace, inline | Helmholtz, inline |
|---|---|---|---|---|---|---|
| 1 | 8.45 | 3.69 | 3.41 | 5.60 | 5.51 | 5.30 |
| 2 | 8.41 | 3.74 | 3.40 | 8.72 | 8.10 | 7.98 |
| 3 | 8.44 | 3.81 | 3.42 | 10.03 | 7.81 | 7.01 |
| 4 | 8.54 | 3.76 | 3.37 | 10.67 | 7.20 | 6.09 |
| 5 | 8.55 | 3.76 | 3.27 | 10.80 | 6.81 | 5.90 |
| 6 | 8.46 | 3.78 | 3.32 | 10.86 | 6.22 | 5.28 |
| 7 | 8.59 | 3.59 | 3.13 | 10.78 | 5.79 | 4.30 |

Throughput in GDoF/s on 64 CPU cores.

| Q | Mass, pre | Laplace, pre | Helmholtz, pre | Mass, inline | Laplace, inline | Helmholtz, inline |
|---|---|---|---|---|---|---|
| 1 | 8.77 | 3.86 | 3.44 | 5.71 | 5.46 | 5.52 |
| 2 | 8.79 | 3.83 | 3.43 | 8.36 | 6.09 | 5.68 |
| 3 | 8.75 | 3.90 | 3.46 | 9.88 | 4.70 | 4.35 |
| 4 | 8.78 | 3.84 | 3.44 | 10.44 | 4.63 | 3.98 |
| 5 | 8.77 | 3.96 | 3.45 | 10.69 | 4.34 | 3.76 |
| 6 | 8.78 | 3.83 | 3.36 | 10.50 | 3.89 | 3.40 |
| 7 | 8.90 | 3.46 | 2.96 | 10.38 | 3.75 | 2.89 |

## Interpretation

### Bytes and flops per dof

Each dof moves the bytes and costs the flops below. Bytes count each argument of
the action once, from the shapes of its arrays. Flops count the additions,
subtractions, multiplications and divisions in each linalg op of the compiled
form, times the size of its iteration space. On the CPU, writing the result also
reads it into cache, which adds 8 bytes per dof.

| Q | Mass, pre | Laplace, pre | Helmholtz, pre | Mass, inline | Laplace, inline | Helmholtz, inline |
|---|---|---|---|---|---|---|
| 1 | 24 B, 26 | 64 B, 92 | 72 B, 119 | 40 B, 247 | 40 B, 376 | 40 B, 403 |
| 2 | 24 B, 38 | 64 B, 128 | 72 B, 167 | 23 B, 259 | 23 B, 412 | 23 B, 451 |
| 3 | 24 B, 50 | 64 B, 164 | 72 B, 215 | 19 B, 271 | 19 B, 448 | 19 B, 499 |
| 4 | 24 B, 62 | 64 B, 200 | 72 B, 263 | 18 B, 283 | 18 B, 484 | 18 B, 547 |
| 5 | 24 B, 74 | 64 B, 236 | 72 B, 311 | 17 B, 295 | 17 B, 520 | 17 B, 595 |
| 6 | 24 B, 86 | 64 B, 272 | 72 B, 359 | 17 B, 307 | 17 B, 556 | 17 B, 643 |
| 7 | 24 B, 98 | 64 B, 308 | 72 B, 407 | 16 B, 319 | 16 B, 592 | 16 B, 691 |

The flops for inline geometry are an upper bound. On the GPU, inline mass runs at
more than the FP64 peak by this count, so the compiled kernels do less work than
the IR describes. Most of it is the Jacobian of each point, and much of that is
likely shared between the points of a pencil once it is scheduled.

### Efficiency

The percentage of the bandwidth reached on the GPU, and for inline geometry also
the percentage of the FP64 peak by the flop count above.

| Q | Mass, pre | Laplace, pre | Helmholtz, pre | Mass, inline | Laplace, inline | Helmholtz, inline |
|---|---|---|---|---|---|---|
| 1 | 100 | 93 | 93 | 91, 63 | 71, 76 | 60, 68 |
| 2 | 94 | 87 | 87 | 68, 86 | 33, 67 | 30, 65 |
| 3 | 95 | 90 | 82 | 64, 102 | 30, 81 | 23, 68 |
| 4 | 98 | 90 | 84 | 62, 114 | 27, 85 | 24, 85 |
| 5 | 99 | 89 | 79 | 45, 89 | 16, 56 | 15, 58 |
| 6 | 97 | 83 | 79 | 55, 115 | 20, 75 | 17, 76 |
| 7 | 98 | 87 | 82 | 64, 140 | 24, 97 | 19, 91 |

The same on 128 CPU cores, with the peak at the base clock.

| Q | Mass, pre | Laplace, pre | Helmholtz, pre | Mass, inline | Laplace, inline | Helmholtz, inline |
|---|---|---|---|---|---|---|
| 1 | 95 | 94 | 96 | 95, 30 | 93, 45 | 90, 46 |
| 2 | 95 | 95 | 96 | 96, 49 | 89, 72 | 87, 78 |
| 3 | 95 | 97 | 96 | 95, 59 | 74, 76 | 67, 76 |
| 4 | 96 | 95 | 95 | 96, 66 | 65, 76 | 55, 72 |
| 5 | 96 | 95 | 92 | 95, 69 | 60, 77 | 52, 76 |
| 6 | 95 | 96 | 94 | 94, 72 | 54, 75 | 46, 74 |
| 7 | 97 | 91 | 88 | 93, 75 | 50, 74 | 37, 64 |

### Precomputed geometry is limited by memory bandwidth

- On the GPU, mass reaches 94–100% of the copy bandwidth, and Laplace and Helmholtz
  reach 79–93%. Throughput hardly changes with degree, because the bytes per dof do
  not.
- On the CPU, every form reaches 88–97% of the triad bandwidth on 128 cores, and
  83–100% on 64. Half the cores are enough to saturate memory.
- The geometry data is most of the traffic. It is 6 or 7 doubles per point for
  Laplace and Helmholtz, against one value read and one written per dof. That is
  why Laplace runs at about a third of the speed of mass.
- So faster code cannot speed these kernels up. Only fewer bytes per dof can.

### Inline geometry trades bytes for flops

Inline geometry reads 24 doubles per cell instead of 1 to 7 per point. That falls
to 16–23 bytes per dof from degree 2. The kernels then do 2–10 times the flops of
precomputed geometry.

- On the CPU, inline mass stays limited by bandwidth at 93–96% and gains 20–26%
  over precomputed from degree 3. Laplace and Helmholtz about double their
  throughput at degrees 2 and 3, to 7–8 GDoF/s. From degree 3 they are limited by
  arithmetic instead. They reach about three quarters of the base-clock peak on 128
  cores. On 64 cores they reach close to all of it, and inline mass goes over it,
  so the cores run above their base clock and the flop count overstates the work
  here too. Doubling the cores from 64 to 128 speeds Laplace at degree 7 up by 1.5
  times, which also points to arithmetic.
- On the GPU, inline geometry matches or slightly beats precomputed for Laplace and
  Helmholtz, except at degrees 5 and 6, and is slower for mass. Laplace and
  Helmholtz gain most at degree 1, where precomputed geometry moves the most bytes.
  Above that the inline kernels are limited by FP64 throughput, so the V100's 8.9
  flops per byte of bandwidth leaves less room to trade bytes for flops than the
  CPU's 16.
- Degree 5 is markedly slower with inline geometry on the GPU, by 25–40% against
  degrees 4 and 6. Precomputed geometry does not show this, so it is likely a
  scheduling effect of the inline kernels at this size. It is not yet explained.

### GPU against CPU

With precomputed geometry the V100 is 3.7–3.9 times faster than 128 CPU cores for
mass and 2.6–3.1 times faster for Laplace and Helmholtz. That is close to the ratio
of their bandwidths, 2.8–3.7 depending on how write-allocate is counted. With
inline geometry the gap narrows to 1.1–2.2 times for Laplace and Helmholtz,
because the CPU has more arithmetic per byte of bandwidth.

## Gather and scatter on the GPU

The same kernels also act on assembled vectors, measured on the GPU on 2026-10-08.
Each kernel gathers its cell values of `x` through the dofmap. The results are then
summed into `y` in one of two ways.

- `CSR`: the kernel writes its values on each cell, and a separate CUDA kernel in
  the driver sums them into `y` through the transposed dofmap in CSR form. It is
  deterministic.
- `atomic`: the kernel adds its values into `y` with atomic adds, after the driver
  zeroes `y`.

All indices are int32, as in DOLFINx. Throughput here counts global dofs, of which
there are fewer than E-vector dofs, by a factor of ((k + 1) / k)³. That is 8 at
degree 1 and 1.5 at degree 7. The percentage is the bandwidth reached against the
798 GB/s of the copy kernel. It counts every array each step reads or writes once,
including the dofmap, the transposed dofmap and the values on cells of `CSR`.

Throughput in GDoF/s, with the percentage of bandwidth in brackets, for
precomputed geometry.

| Q | Mass, CSR | Mass, atomic | Laplace, CSR | Laplace, atomic | Helmholtz, CSR | Helmholtz, atomic |
|---|---|---|---|---|---|---|
| 1 | 1.96 (67) | 4.45 (71) | 1.05 (77) | 1.44 (80) | 0.95 (78) | 1.27 (81) |
| 2 | 3.65 (58) | 7.60 (69) | 2.07 (68) | 3.14 (81) | 1.79 (65) | 2.68 (78) |
| 3 | 6.67 (79) | 9.31 (70) | 3.45 (82) | 3.98 (77) | 2.96 (77) | 3.44 (75) |
| 4 | 7.67 (79) | 10.69 (74) | 4.04 (81) | 4.89 (81) | 3.54 (77) | 4.15 (77) |
| 5 | 8.79 (82) | 11.03 (73) | 4.48 (80) | 4.96 (75) | 3.73 (73) | 3.99 (67) |
| 6 | 8.85 (78) | 11.39 (73) | 4.75 (79) | 5.35 (76) | 4.17 (76) | 4.52 (72) |
| 7 | 9.79 (83) | 11.82 (74) | 5.10 (81) | 5.53 (76) | 4.45 (77) | 4.64 (70) |

The same for inline geometry.

| Q | Mass, CSR | Mass, atomic | Laplace, CSR | Laplace, atomic | Helmholtz, CSR | Helmholtz, atomic |
|---|---|---|---|---|---|---|
| 1 | 1.37 (69) | 2.27 (72) | 1.15 (58) | 1.73 (55) | 1.05 (53) | 1.47 (47) |
| 2 | 3.10 (48) | 6.03 (52) | 2.15 (33) | 3.20 (28) | 2.00 (31) | 2.80 (24) |
| 3 | 6.16 (64) | 9.01 (55) | 3.96 (41) | 4.92 (30) | 3.22 (34) | 3.80 (23) |
| 4 | 7.12 (62) | 10.84 (58) | 4.53 (39) | 5.67 (30) | 4.12 (36) | 5.11 (27) |
| 5 | 7.02 (55) | 9.80 (49) | 3.44 (27) | 3.96 (20) | 3.04 (24) | 3.49 (18) |
| 6 | 7.88 (58) | 10.73 (53) | 4.29 (32) | 4.99 (24) | 3.86 (28) | 4.47 (22) |
| 7 | 9.57 (67) | 12.49 (60) | 5.28 (37) | 6.30 (30) | 4.15 (29) | 4.69 (23) |

- With precomputed geometry both scatters reach 58–83% of the bandwidth, so they
  are limited by memory like the action on E-vectors.
- Atomic adds are faster than the CSR scatter in every case, by 1.04–2.3 times. The
  CSR scatter writes and reads back 16 bytes per value on cells, which costs most
  at low degree, where each global dof has the most values on cells. At degrees 1
  and 2 the CSR scatter takes 55–61% of the time for mass with precomputed geometry,
  and 45% at degree 7.
- Atomic adds cost 24 bytes of `y` per global dof, for zeroing it then reading and
  writing it. Writing `y` once would cost 8.
- Throughput rises with degree, unlike on E-vectors, because the ratio of values on
  cells to global dofs falls. Mass at degree 7 with atomic adds runs at 11.8 GDoF/s,
  which is 17.6 G values on cells per second. On E-vectors it runs at 32.7. The
  difference is the dofmap, the reads of `x` through it and the traffic of `y`.
- Inline geometry beats precomputed for Laplace and Helmholtz with atomic adds at
  every degree except 5 and 6, by up to 23% for Helmholtz at degree 4. The extra
  traffic of the scatter makes the bytes it saves worth more. It reaches only
  18–30% of the bandwidth from degree 2, so these kernels are limited by
  arithmetic, as on E-vectors. Degree 5 is slow again.

## Reproducing

From a shell with the Spack environment and the virtual environment active.

```sh
cmake -G Ninja -S bench/gpu -B build/bench_gpu \
  -DFFCY_FORMS="mass;laplace;helmholtz" -DFFCY_DEGREES="1;2;3;4;5;6;7" \
  -DFFCY_GEOMETRIES="precomputed;inline" -DFFCY_CELLS_PER_BLOCK="1;2;4;8;16"
cmake --build build/bench_gpu -- -k 0
python -m bench.sweep build/bench_gpu --seconds 1 --num-dofs 3.2e7

cmake -DFFCY_VECTORS="evector;assembled;fused" build/bench_gpu
cmake --build build/bench_gpu -- -k 0
python -m bench.sweep build/bench_gpu --filter _assembled --seconds 1 --num-dofs 3.2e7
python -m bench.sweep build/bench_gpu --filter _fused --seconds 1 --num-dofs 3.2e7

cmake -G Ninja -S bench/cpu -B build/bench_cpu \
  -DFFCY_FORMS="mass;laplace;helmholtz" -DFFCY_DEGREES="1;2;3;4;5;6;7" \
  -DFFCY_GEOMETRIES="precomputed;inline" -DFFCY_CELLS_PER_BLOCK="4;8"
cmake --build build/bench_cpu
OMP_NUM_THREADS=128 OMP_PROC_BIND=spread OMP_PLACES=cores \
  python -m bench.sweep build/bench_cpu --seconds 1 --num-dofs 6.4e7
```

The GPU build stops at the kernels that use too much shared memory, so it needs
`-k 0` to build the rest.
