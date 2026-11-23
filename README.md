# Eductor dashboard

Private source replica of the analytical exposure lab. This repository does not change the hosted original.

## Open the simulator

The running simulator is the original site:

https://rev2-eductor-dashboard.a-spector.workers.dev/analytical-exposure

This GitHub repository is the readable source, not a hosted app.

## Layout

- `dashboard/analytical-exposure.html` loads the lab.
- `dashboard/src/model/` is the flow, particle, collision, and calibration code.
- `dashboard/src/ui/lab.js` is the experiment screen.
- `dashboard/src/data/` holds CAD meshes and render images.
- `dashboard/src/vendor/` holds Three.js and the orbit controls.
- `model-assets/exposure-worker.js` runs the particle histories.
- `downloads/` holds the linked reports.

Serve the `dashboard` directory and `model-assets` from the same site root so `/model-assets/exposure-worker.js` resolves.
