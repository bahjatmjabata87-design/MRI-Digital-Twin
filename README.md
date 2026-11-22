# MRI Digital Twin – Safety & Contrast Lab (v1.0)

An interactive educational simulator for radiography and medical students. It runs in any modern browser with no installation.

**Live demo:** open `index.html`, or use the GitHub Pages link in this repository's About section.

## Modules

1. **MRI suite and safety zones (3D)**
   - A 3D MRI suite divided into the four ACR MR safety zones (I–IV), with an explanation of each zone.
   - An illustrative 5-gauss (0.5 mT) line around the magnet.
   - A projectile test: bring objects into Zone IV and see how MR Unsafe (ferromagnetic), MR Conditional and MR Safe items behave.
2. **Contrast lab (TR / TE)**
   - Adjust repetition time (TR), echo time (TE) and field strength (1.5 T or 3 T).
   - A simplified axial brain model (fat, CSF, grey and white matter) updates live, with T1-weighted, T2-weighted and proton-density presets.
   - Uses the basic spin-echo signal equation: `S = PD × (1 − e^(−TR/T1)) × e^(−TE/T2)`.

## Built with

- HTML, CSS and JavaScript in one file
- [three.js](https://threejs.org/) r128 for the 3D scene

## Scientific basis

- Zone definitions: American College of Radiology (ACR) Manual on MR Safety, 2024.
- Object labels: ASTM F2503 terms (MR Safe, MR Conditional, MR Unsafe).
- Brain relaxation times at 3 T: Wansapura JP et al. J Magn Reson Imaging. 1999;9(4):531–538.
- Fat at 3 T: Radiat Oncol. 2013;8:51, Table 1.
- Other values (1.5 T brain, CSF, fat at 1.5 T, proton densities) are rounded approximations of published values. Values vary with measurement method, region and age.

## Disclaimer

For education only. This is not a medical device, does not model any specific scanner, and must not be used for clinical decisions, patient screening or MR site planning. Always follow your site's MR safety policy and the labelling of each item.

## Roadmap

- v1.1: gradient-echo and inversion-recovery (FLAIR, STIR) contrast
- v1.2: k-space and image formation module
- v1.3: MR safety screening quiz

## Author

Bahjat Jabatah, MPH – MRI/CT Technologist and Health Researcher. Feedback and suggestions are welcome through Issues.

## License

MIT License – see `LICENSE`.
