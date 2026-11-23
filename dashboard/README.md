# Eductor dashboard

Private local replica of the analytical exposure lab. The original hosted dashboard is not modified by this repository.

## Layout

- `analytical-exposure.html` loads the lab in the browser.
- `src/model/` is the readable flow, particle, collision, and calibration code.
- `src/ui/lab.js` is the experiment screen.
- `src/data/` holds CAD meshes and render images.
- `src/vendor/` holds Three.js and the orbit controls.
- `../model-assets/exposure-worker.js` runs the particle histories.

Open `analytical-exposure.html` from a local web server so the worker and scripts can load.
