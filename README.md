# Eductor dashboard

Private local replica of the analytical exposure lab. The original hosted dashboard is not modified by this repository.

## Layout

- `dashboard/analytical-exposure.html` loads the lab in the browser.
- `dashboard/src/model/` is the readable flow, particle, collision, and calibration code.
- `dashboard/src/ui/lab.js` is the experiment screen.
- `dashboard/src/data/` holds CAD meshes and render images.
- `dashboard/src/vendor/` holds Three.js and the orbit controls.
- `model-assets/exposure-worker.js` runs the particle histories.

Open the HTML from a local web server so the worker and scripts can load.
