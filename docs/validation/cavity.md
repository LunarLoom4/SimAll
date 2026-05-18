# Validation — Lid-Driven Cavity

3-D lid-driven cavity, $1 \times 1 \times 1$ m, top lid at
$U = 1$ m/s, all other walls no-slip.

| Re | Mesh | Scheme | Reference | Tolerance |
|---|---|---|---|---|
| 100  | $32^3$ uniform hex | SIMPLE, UD | Ghia 1982 | u-mid centreline within 1 % |
| 400  | $64^3$             | SIMPLE, MUSCL | Ghia 1982 | 1 % |
| 1000 | $96^3$             | SIMPLE, MUSCL | Ghia 1982 | 1 % |
| 3200 | $128^3$            | Coupled, MUSCL | Ghia 1982 | 2 % |

Run from `applications/simall_cavity`. Reference data live in
`data/fields/ghia1982.csv`.
