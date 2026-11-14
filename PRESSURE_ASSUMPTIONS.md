# Values used to turn Cp into pressure

This records every number that went into the pressures in `out/e423_points.csv` and
`out/e423_buckets.csv` (E423, α = 12°), and says which were given, which were
defaults, and which were assumed.

## Short answers

- **Sea-level constants: yes.** Density is the standard-atmosphere sea-level value,
  ρ = 1.225 kg/m³. It is the tool's default and was never overridden.
- **Characteristic length: none was used.** Converting Cp to pressure needs only the
  dynamic pressure, so no chord enters the pressure numbers. The reference run was
  not given a chord at all.
- **Airspeed came from the Mach number, not the Reynolds number.** V = 34 m/s is
  Mach 0.1 at the sea-level speed of sound (340.3 m/s). Nothing was worked backwards
  from Re to get the pressures.
- **One thing to fix:** Re = 9×10⁵, V = 34 m/s and your 25.5 in chord cannot all be
  true at sea level. See [the inconsistency](#the-inconsistency-to-resolve) below.

## The conversion

```
p = n · K · q · Cp          q = ½ ρ V²
```

Command that produced the outputs:

```bash
python3 cp_visualizer.py e423.dat e423_a12.cp --V 34 --spacing 3 --thickness 0.125 --buckets 6 --bucket-spacing adaptive --alpha 12 --save out/e423
```

| Quantity | Value | Where it came from |
|---|---|---|
| Air density ρ | 1.225 kg/m³ | tool default, sea level at 15 °C; not passed on the command line |
| Airspeed V | 34 m/s (76.1 mph) | `--V 34`; chosen as Mach 0.1 at sea level, **not a confirmed flight speed** |
| Dynamic pressure q | ½ · 1.225 · 34² = **708.05 Pa** (0.1027 psi) | computed |
| Rib spacing | 3 in | `--spacing 3`, your value |
| Rib thickness | 0.125 in | `--thickness 0.125`, your value |
| Edge factor K | 3 / 0.125 = **24** | computed; interior rib (an end rib would be 12) |
| Load factor n | **1** | tool default; no g-load or safety factor applied |
| Multiplier n·K·q | **16 993.2 Pa per unit Cp** (2.4647 psi) | computed |
| Pa to psi | 6894.757 Pa per psi | constant in the tool |
| Chord | not used | `--chord` was not passed |
| Mach correction | none | `--mach` not passed, and the formula has no Mach term |

Check against the output: the first point has Cp = −0.05684, and
16 993.2 × (−0.05684) = −965.89 Pa, which is the value in `out/e423_points.csv`.
The suction peak is Cp = −3.2132, giving −54 603 Pa (−7.92 psi) on the rib edge.

## What the pressure numbers mean

- **They are rib-edge pressures, not air pressures.** K = 24 concentrates the load
  from a 3 in strip of covering onto a 1/8 in rib edge. The pressure the air actually
  exerts on the skin is q·Cp, 24 times smaller: 708.05 Pa per unit Cp, so about
  −2275 Pa (−0.33 psi) at the suction peak.
- **They are relative to freestream static pressure.** Atmospheric pressure
  (101 325 Pa at sea level) is not added, because it acts all the way round a closed
  section and cancels.
- **They are for 1 g level conditions at this α.** No manoeuvre load factor or
  safety factor is included.

## Values that sit behind the Cp file

These were XFOIL inputs. They shape the Cp distribution but are not used again in the
conversion.

| Quantity | Value | Note |
|---|---|---|
| Reynolds number | 9×10⁵ | viscous run; sets the boundary layer and therefore the Cp shape |
| Mach number | 0.1 | XFOIL applies its own compressibility correction to Cp, about +0.5 % at this Mach |
| Angle of attack | 12° | attached flow, trailing-edge Cp ≈ −0.06 |
| Chord in XFOIL | 1 (normalized) | coordinates are x/c, y/c |

Sea-level values implied by "Mach 0.1 ≈ 34 m/s" but never typed into the tool:
temperature 15 °C (288.15 K), speed of sound 340.3 m/s, dynamic viscosity
1.789×10⁻⁵ Pa·s, kinematic viscosity ν = 1.461×10⁻⁵ m²/s.

## The inconsistency to resolve

Reynolds number ties speed and chord together: Re = V·c / ν.

| Take as fixed | Then | Result |
|---|---|---|
| Re = 9×10⁵ and V = 34 m/s | chord = Re·ν / V | 0.387 m = **15.2 in** |
| V = 34 m/s and chord = 25.5 in | Re = V·c / ν | **1.51×10⁶** |
| Re = 9×10⁵ and chord = 25.5 in | V = Re·ν / c | **20.3 m/s** (45 mph), Mach 0.060 |

The 15.2 in figure matches the `--chord 0.386` that appears in the tool's own example
command and in the beginner's guide, so the original examples were consistent with a
15.2 in chord. The discretizer now uses your 25.5 in chord, which does not match.

This matters most in the third case. Pressure scales with V², so at 20.3 m/s the
dynamic pressure is 252.3 Pa instead of 708.05 Pa, and **every pressure in the
outputs would be 0.356 times its current value**, about 2.8 times smaller.

What decides it is your real flight speed at the design condition:

- If the wing really flies at about 34 m/s with a 25.5 in chord, the pressures stand
  as they are, but XFOIL should be re-run at Re ≈ 1.5×10⁶ so the Cp shape matches.
  The change in Cp from 9×10⁵ to 1.5×10⁶ is usually small for attached flow, but it
  moves stall.
- If Re = 9×10⁵ is the right number for a 25.5 in chord, re-run the tool with
  `--V 20.3` and the pressures drop by the factor above.

## Still open

- Real flight speed at the design case (34 m/s is a placeholder from Mach 0.1).
- Load factor and safety factor (currently n = 1; typical RC design is 3 to 5 g).
- Flying altitude and temperature, if not near sea level at 15 °C. Density changes q
  in direct proportion.
