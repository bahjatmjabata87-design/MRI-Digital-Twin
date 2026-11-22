# MRI Digital Twin – Safety & Contrast Lab (v1.0)

An interactive educational simulator for radiography and medical students. It runs in any modern browser with no installation.

**Live demo:** https://bahjatmjabata87-design.github.io/MRI-Digital-Twin/

An 80-second video tour (`demo.mp4`) plays at the top of the page.

## Modules

1. **Contrast lab (TR / TE)**
   - Adjust repetition time (TR), echo time (TE) and field strength (1.5 T or 3 T), with T1-weighted, proton-density and T2-weighted presets.
   - Optional realistic mode: brightness follows the true signal and Rician noise is added, so weak settings look dark and grainy and 3 T looks cleaner (simplified SNR model).
   - A brain MRI image is computed live from real brain anatomy: the grey matter, white matter and CSF maps of the MNI ICBM152 2009a template, across four axial slices.
   - Uses the basic spin-echo signal equation: `S = PD × (1 − e^(−TR/T1)) × e^(−TE/T2)`.
2. **MRI suite and safety zones (3D)**
   - A 3D MRI suite divided into the four ACR MR safety zones (I–IV), with an explanation of each zone.
   - An illustrative 5-gauss (0.5 mT) line around the magnet.
   - A projectile test with 3D models (steel and aluminium oxygen cylinders, scissors, wheelchair, infusion pump): bring them into Zone IV and see how MR Unsafe (ferromagnetic), MR Conditional and MR Safe items behave.

## Built with

- HTML, CSS and JavaScript in one file
- [three.js](https://threejs.org/) r128 for the 3D scene

## Scientific basis

- Zone definitions: American College of Radiology (ACR) Manual on MR Safety, 2024.
- Object labels: ASTM F2503 terms (MR Safe, MR Conditional, MR Unsafe).
- Brain relaxation times at 3 T: Wansapura JP et al. J Magn Reson Imaging. 1999;9(4):531–538.
- Brain anatomy: MNI ICBM152 2009a nonlinear symmetric template (Fonov VS et al. NeuroImage 2011;54(1):313–327; NeuroImage 2009;47(Suppl 1):S102), grey and white matter maps as distributed with nilearn. CSF is estimated as brain volume not classed as grey or white matter. See `NOTICE`.
- Other values (1.5 T brain, CSF, proton densities) are rounded approximations of published values. Values vary with measurement method, region and age.

The images are an average of many adult scans, not a single patient.

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
