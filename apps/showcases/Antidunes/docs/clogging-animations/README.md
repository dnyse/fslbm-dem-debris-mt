# Clogging & equilibrium demo animations

MP4 animations referenced from the project's top-level README and from
`defense_presi` (as Theora `.avi` copies, see that repo's `videos/README.md`).
All 1 frame per simulated snapshot, `-r 30` output (each source frame held
for 1 second of playback). `*-compact-*` are a cropped/scaled layout for
embedding (web pages, slides); the plain `*-1fps.mp4` files (where present)
are the full-resolution multi-view renders.

| file | run | source (final_anim/) |
|------|-----|------------------------|
| `eq-sphere-{,compact-}1fps.mp4`   | sphere equilibrium, ρ=0.25/0.50/0.75 side by side | `eq_runs/spheres/` (strict-Bogner-classification runs 4408737/4408739/4408741) |
| `eq-cylinder-{,compact-}1fps.mp4` | cylinder equilibrium, ρ=0.25/0.50/0.75 side by side | `eq_runs/cyl/` (strict-Bogner-classification runs 4408736/4408738/4408740) |
| `eq-cuboid-{,compact-}1fps.mp4`   | cuboid equilibrium, ρ=0.25/0.50/0.75 side by side | `eq_runs/cuboid/` (water-viscosity reruns 4351201/4351202/4351203) |
| `clog-406rods-{,compact-}1fps.mp4` | 406 rods jam the hopper gap | `4293495_hopper_406rods_startblock_20261004/` |
| `clog-406boxes-compact-1fps.mp4`   | 406 cuboids in the hopper gap | `4297652_hopper_406boxes_startblock_20261004/` |
| `clog-mixed-compact-1fps.mp4`      | mixture of shapes jams the hopper gap | `4293831_hopper_mixed_startblock_20261004/` |
| `clog-twopillar-compact-1fps.mp4`  | 812 rods at the two-pillar layout (porous, non-jamming) | `4300948_twopillar_812rods_startblock_20261005/` |

The clogging clips are each 3 camera views (bird/front/side, or front/side/top
depending on the run) composited side by side; the equilibrium clips are the
3 density ratios side by side, left to right.
