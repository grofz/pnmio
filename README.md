# PNM Fortran Library

A small Fortran library for reading and writing PNM image files (PPM, PGM, PBM) and basic colormap handling.

## Features

- **readpnm** – Read PNM images (PPM, PGM, PBM).
- **writeppm** – Write PPM images (color images).
- **writepgm** – Write PGM (grayscale) or PBM (bitmap) images.
- **colormap** – Get RGB colormap arrays.
- **assign\_colormap** – Map scalar fields to RGB images using a colormap.

## Use examples

Read a PPM image

```fortran
  use pnmio_module
  :
  integer, allocatable :: aa(:,:,:)
  call readpnm('file.ppm', aa)
```

Read a PGM or PBM image

```fortran
  use pnmio_module
  :
  integer, allocatable :: aa(:,:)
  call readpnm('file.pgm', aa)
```

Write PPM image (in plain format)

```fortran
  use pnmio_module
  :
  integer, dimension(300, 100) :: rr, gg, bb
  :
  rr = ...
  gg = ...
  bb = ...
  call writeppm('output.ppm', rr, gg, bb, is_plain=.true.)
```

Write PBM image

```fortran
  use pnmio_module
  :
  integer :: aa(300, 100)
  :
  aa = ... ! set 0 or 1
  call writepgm('bitmap.pbm', aa, mx=1)
```

---

Visualize a 2D field by making a PPM image

```Fortran
  use pnmio_module
  use iso_fortran_env, only : dp=>real64
  :
  integer, parameter :: NX = 320, NY = 240
  real(dp) :: uu(NX, NY)
  integer, dimension(NX, NY) :: rr, gg, bb

  call random_number(uu)
  call assign_colormap(uu, rr, gg, bb, 30, idcolormap=CM_TURBO)
  call writeppm('randomized.ppm', rr, gg, bb)
```

## Full working examples

Iteratively solve PDE for heat conduction

```sh
$ gfortran -o a.out pnmio.f90 example/main_heat.f90
$ ./a.out
$ okular test0.ppm test.ppm
```

Visualize the colormaps

```sh
$ gfortran -o b.out pnmio.f90 example/main_showcolormap.f90
$ ./b.out
$ okular test_colormap*.ppm
```

---

## Interfaces and Usage

### 1. `readpnm`

```fortran
call readpnm(filename, aa [, ierr])
```

**Inputs:**

- `FILENAME` (character) – File name of the PNM image to read.

**Outputs:**

- `AA` (integer, allocatable) –
  - Rank-3 array for color images (PPM), with dimensions:\
    `AA(1:3, width, height)` → Red, Green, Blue components.
  - Rank-3 array for grayscale/bitmap images (PGM/PBM), with first index set to 1:\
    `AA(1, width, height)`
  - Alternatively, a rank-2 array can be used for grayscale/bitmap images:\
    `AA(width, height)`
- `IERR` (integer, optional) – Non-zero indicates an error.

### 2. `writeppm`

```fortran
call writeppm(filename, rr, gg, bb [, mx] [, is_plain] [, ierr])
call writeppm(filename, aa [, mx] [, is_plain] [, ierr])
```

**Inputs:**

- `FILENAME` (character) – Output file name.
- `RR`, `GG`, `BB` (integer, rank-2 arrays) – Color components.
- `AA` (integer, rank-3 array) – Color components where:
  - AA(1,:,:) – Red
  - AA(2,:,:) – Green
  - AA(3,:,:) – Blue
- `MX` (integer, optional) – Maximum color value (default: auto-set to 255 or 65535).
- `IS_PLAIN` (logical, optional) – If `.TRUE.`, writes ASCII (plain) format. Default is binary (raw).

**Outputs:**

- `IERR` (integer, optional) – Non-zero indicates an error.

### 3. `writepgm`

```fortran
call writepgm(filename, aa [, mx] [, is_plain] [, ierr])
```

**Inputs:**

- `FILENAME` (character) – Output file name.
- `AA` (integer, rank-2 array) – Image data.
- `MX` (integer, optional) –
  - `MX = 1` -> Write PBM file (bitmap, black & white).
  - `MX not equal to 1` -> Write PGM file (grayscale). Default: auto-set to 255 or 65535.
- `IS_PLAIN` (logical, optional) – If `.TRUE.`, writes ASCII (plain) format. Default is binary (raw).

**Outputs:**

- `IERR` (integer, optional) – Non-zero indicates an error.

---

### 4. `colormap`

```fortran
call colormap(red, green, blue [, idcolormap])
```


**Inputs:**

- `IDCOLORMAP` (integer, optional) – Selects the colormap:

| ID | Colormap                 |
| -- | ------------------------ |
| 1  | CM\_RAINBOW (Jet-style)  |
| 2  | CM\_VIRIDIS (Default)    |
| 3  | CM\_TURBO (Improved Jet) |

**Outputs:**

- `RED`, `GREEN`, `BLUE` (integer, rank-1 arrays) – Colormap RGB components.

- First and last elements of the output arrays are black and white, respectively.
- Intermediate elements correspond to equidistant colors from the selected colormap.

### 5. `assign_colormap`

```fortran
call assign_colormap(aa, rr, gg, bb, n [, amin], [, amax], [, idcolormap])
call assign_colormap(uu, rr, gg, bb, n [, umin], [, umax], [, idcolormap])
```

**Inputs:**

- `AA` (integer, rank-2 array) – Input scalar field (integer).
- `UU` (real, rank-2 array) – Input scalar field (real).
- `N` (integer) – Number of colors in the colormap.
- `AMIN`, `AMAX` (integer, optional) – Min/Max for integer input mapping.
- `UMIN`, `UMAX` (real, optional) – Min/Max for real input mapping.
- `IDCOLORMAP` (integer, optional) – Colormap ID (same as in `COLORMAP`).

**Outputs:**

- `RR`, `GG`, `BB` (integer, rank-2 arrays) – RGB components of the output image (0–255).

**Notes:**

- Values outside the given limits are mapped to black (below) and white (above).
- If limits are not provided, they default to the min/max of the data.

---

## Contact

grofz@vscht.cz (Z. Grof)

## License

This library is provided as-is, with no warranty. Free to use for academic or personal projects.

