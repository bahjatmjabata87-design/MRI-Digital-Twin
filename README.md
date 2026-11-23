# MRI Digital Twin – Contrast, Planning & Safety Lab (v1.1)

An interactive educational simulator for radiography and medical students. It runs in any modern browser with no installation.

**Live demo:** https://bahjatmjabata87-design.github.io/MRI-Digital-Twin/

An 80-second video tour (`demo.mp4`) plays at the top of the page. The planning tab shows its own 103-second tour (`planning-tour.mp4`).

## Modules

1. **Contrast lab (TR / TE)**
   - Adjust repetition time (TR), echo time (TE) and field strength (1.5 T or 3 T), with T1-weighted, proton-density and T2-weighted presets.
   - Optional realistic mode: brightness follows the true signal and Rician noise is added, so weak settings look dark and grainy and 3 T looks cleaner (simplified SNR model).
   - A brain MRI image is computed live from real brain anatomy: the grey matter, white matter and CSF maps of the MNI ICBM152 2009a template, across four axial slices.
   - Uses the basic spin-echo signal equation: `S = PD × (1 − e^(−TR/T1)) × e^(−TE/T2)`.
2. **Positioning and planning (brain) – new in v1.1**
   - Position a virtual patient: head first or feet first, supine or prone, head coil or body coil, foam pads, hearing protection, and the laser landmark (glabella, nasion, vertex, chin or sternal notch). The side view is drawn from the same 3D head, so the laser lands on the chosen landmark. Move the table to the isocentre to acquire three-plane localizers.
   - Plan axial, coronal or sagittal slices on the sagittal, coronal and axial localizers: drag to move the slice block and drag the handle to angle it (AC–PC line guide available). Each new patient has a random head tilt to correct.
   - Set TR, TE, echo train length, NEX, number of slices, thickness, gap, FOV, matrix and phase-encoding direction, at 1.5 T or 3 T. Coverage, voxel size, relative SNR, slices per TR and scan time update live.
   - Scan and scroll through the resulting slices, resliced live from a 3D brain volume. See partial volume, pixelation, noise, wrap-around and motion ghosts.
   - "Check my scan" scores the positioning, plan and protocol against a task (routine axial T2, coronal T2 or sagittal T1). "Show me" demonstrates a correct plan.
3. **MRI suite and safety zones (3D)**
   - A 3D MRI suite divided into the four ACR MR safety zones (I–IV), with an explanation of each zone.
   - An illustrative 5-gauss (0.5 mT) line around the magnet.
   - A projectile test with 3D models (steel and aluminium oxygen cylinders, scissors, wheelchair, infusion pump): bring them into Zone IV and see how MR Unsafe (ferromagnetic), MR Conditional and MR Safe items behave.

## Built with

- HTML, CSS and JavaScript (`index.html`, `planning.js`); 3D head data in `brain3d.js` (4-bit grey matter, white matter, CSF, fat, soft tissue and eye maps on a 2 mm grid)
- [three.js](https://threejs.org/) r128 for the 3D scene

## Scientific basis

- Zone definitions: American College of Radiology (ACR) Manual on MR Safety, 2024.
- Object labels: ASTM F2503 terms (MR Safe, MR Conditional, MR Unsafe).
- Brain relaxation times at 3 T: Wansapura JP et al. J Magn Reson Imaging. 1999;9(4):531–538.
- Brain anatomy: MNI ICBM152 2009a nonlinear symmetric template (Fonov VS et al. NeuroImage 2011;54(1):313–327; NeuroImage 2009;47(Suppl 1):S102), grey and white matter maps as distributed with nilearn. CSF is estimated as brain volume not classed as grey or white matter. See `NOTICE`.
- Planning rules: axial slices parallel to the AC–PC line from vertex to foramen magnum; coronal perpendicular to AC–PC; sagittal parallel to the midline. Scan time = TR × phase steps × NEX ÷ ETL; SNR ∝ voxel volume × √(phase steps × NEX). AC and PC were located approximately on the template.
- Scalp fat at 3 T: de Bazelaire CMJ et al. Radiology. 2004;230(3):652–659.
- Other values (1.5 T brain, CSF, proton densities) are rounded approximations of published values. Values vary with measurement method, region and age.

The skull, scalp, face, eyes, neck and spinal cord in the planning lab are simplified illustrative shapes around the template brain. The images are an average of many adult scans, not a single patient.

## Disclaimer

For education only. This is not a medical device, does not model any specific scanner, and must not be used for clinical decisions, patient screening or MR site planning. Always follow your site's MR safety policy and the labelling of each item.

## Roadmap

- v1.2: gradient-echo and inversion-recovery (FLAIR, STIR) contrast
- v1.3: k-space and image formation module
- v1.4: more planning anatomy (spine, knee) and an MR safety screening quiz

## Author

Bahjat Jabatah, MPH – MRI/CT Technologist and Health Researcher. Feedback and suggestions are welcome through Issues.

## License

MIT License – see `LICENSE`.
