# Planetary Alignment

[Planetary Alignment Visualization](https://planetary-align.streamlit.app/)

A computational tool for identifying multi-body celestial alignments using astronomical coordinate systems and parallel processing.

## Table of Contents
- [Key Features](#key-features)
- [Astronomical Concepts Implemented](#astronomical-concepts-implemented)
- [Code Structure & Algorithms](#code-structure--algorithms)
- [Customization Options](#customization-options)
- [Interpreting Results](#interpreting-results)
- [Limitations & Approximations](#limitations--approximations)
- [References](#references)

## Key Features

### 1. Multi-Mode Alignment Detection
- **Pairwise Angular Separation**: Traditional approach checking maximum great-circle distance between any two bodies in ecliptic coordinates
- **Longitude Clustering**: Alternative method considering only ecliptic longitude spread (ignores latitude differences)

### 2. Observability Constraints
- Local altitude threshold filtering (horizon clearance)
- Real-time altazimuth coordinate transformation
- Visibility window optimization

### 3. High Temporal Resolution
- Configurable time steps (1-24 hours)
- Parallelized date processing for efficiency
- Five-year search horizon

## Astronomical Concepts Implemented

### 1. Coordinate Systems
- **Geocentric True Ecliptic**: 
  - Fundamental plane: Earth's orbital plane (ecliptic)
  - Coordinates: Celestial longitude (λ) and latitude (β)
  - Uses VSOP87 theory via JPL DE432s ephemeris
- **Topocentric Horizontal**:
  - Observer-centric coordinates (altitude/azimuth)
  - Accounts for Earth rotation and observer location
  - Uses Astropy's ITRS transformation chain

### 2. Angular Separation Mathematics
- Great-circle distance formula:
$$
  \Delta\sigma = \arccos(\sin\beta_1\sin\beta_2 + \cos\beta_1\cos\beta_2\cos\Delta\lambda)
$$
- Ecliptic longitude spread calculation:
$$
  \text{Range} = \max(\lambda_i) - \min(\lambda_i) \quad \text{(mod 360°)}
$$

### 3. Planetary Motion Considerations
- Inner planet orbital periods:
  - Mercury: 88 days
  - Venus: 225 days
  - Moon: 27.3 days (sidereal)
- Outer planet synodic periods:
  - Jupiter: 398 days
  - Saturn: 378 days

## Code Structure & Algorithms

### Core Functions

1. **`get_ecliptic_positions()`**
   - Ephemeris: JPL DE432s (1 arcsec accuracy)
   - Transformations: ICRS → GCRS → GeocentricTrueEcliptic
   - Output: Dictionary of λ/β positions

2. **`check_ecliptic_alignment()`**
   - Method dispatch for pairwise/longitude checks
   - Circular statistics for longitude wrapping

3. **`get_horizontal_positions()`**
   - Location-specific transformations:
     - EarthLocation → ITRS → AltAz
   - Atmospheric refraction: Not modeled

4. **Parallel Processing**
   - ThreadPoolExecutor for date iteration
   - Early termination on first valid alignment

### Critical Algorithms

```python
def check_alignment_for_date(...):
    # Astronomical pipeline
    positions = get_ecliptic_positions(bodies, date)
    aligned = check_method(positions)
    if aligned and check_visibility:
        altaz = get_horizontal_positions(...)
        visible = all(alt > threshold)
    return (date, aligned & visible)
```

## Customization Options

### 1. Alignment Parameters
| Parameter            | Range      | Astronomical Rationale                           |
| -------------------- | ---------- | ------------------------------------------------ |
| Separation Threshold | 5°-45°     | Moon's angular diameter ≈ 0.5°, conjunction < 5° |
| Time Step            | 1-24 hours | Mercury's daily motion ≈ 4°, Moon ≈ 13°          |
| Minimum Altitude     | -20°-30°   | Atmospheric extinction limits (0° horizon)       |

### 2. Detection Tradeoffs
- **Conservative Settings** (5°, 1h steps):
  - Likely miss events but high precision
- **Permissive Settings** (30°, 6h steps):
  - More false positives but better coverage

## Interpreting Results

### Sky Plot Analysis
- Polar projection (altitude = radial axis)
- Azimuth clockwise from North
- Color-coded body markers

### Physical Interpretation
- **Ecliptic Alignment**: Solar system plane projection
- **Local Sky Position**: Actual visual appearance
- **Historical Context**: Notable alignments:
  - 2000 BC: 5-planet alignment
  - 2020: Jupiter-Saturn great conjunction (0.1°)

## Limitations & Approximations

1. **Geocentric Assumptions**
   - Ignores parallax for Moon (up to 1° error)
   - Assumes instantaneous Earth position

2. **Time Resolution Limits**
   - Mercury's maximum angular speed: 4.1°/day
   - 6h steps → ±0.5° position uncertainty

3. **Atmospheric Effects**
   - No twilight/sun position checks
   - Refraction not modeled below 0°

4. **Ephemeris Limitations**
   - DE432s accuracy (1550-2650 CE)
   - No nongravitational forces for small bodies

## References

1. Astropy Collaboration (2022). "Astropy v5.1" 
   - ICRS/GCRS transformation models
2. Folkner et al. (2014). "JPL Planetary Ephemerides"
   - DE432s specification
3. Meeus (1991). "Astronomical Algorithms"
   - Angular separation mathematics
4. USNO Circular 179 (2010)
   - Coordinate system definitions

---

**License**: MIT Open Source  
**Contact**: astronomy-toolkit@example.com  
**Version**: 1.2 (2024-03-20)
```
