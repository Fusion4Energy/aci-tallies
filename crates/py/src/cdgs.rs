//! D1SUNED CDGS files.

use std::path::PathBuf;
use std::sync::Arc;

use aci_tallies::cdgs::{Cdgs, CdgsInfo};
use pyo3::prelude::*;

use crate::common::{io_error, read_error, warnings, PyAxis, PyReadWarning};
use crate::geometry::{geometry, PyGeometry};
use crate::meshtal::PyMeshTally;

/// Open a CDGS file (read whole; each mesh is parsed when read).
#[pyfunction]
pub fn open_cdgs(py: Python<'_>, path: PathBuf) -> PyResult<PyCdgsFile> {
    let file = py.allow_threads(|| Cdgs::open(&path)).map_err(io_error)?;
    Ok(PyCdgsFile {
        file: Arc::new(file),
    })
}

/// A D1SUNED Common Decay Gamma Source file: a mapping from mesh id to its
/// source as a `MeshTally`.
#[pyclass(name = "CdgsFile", module = "aci_tallies.cdgs", frozen)]
pub struct PyCdgsFile {
    file: Arc<Cdgs>,
}

#[pymethods]
impl PyCdgsFile {
    #[getter]
    fn path(&self) -> Option<PathBuf> {
        self.file.path().map(|p| p.to_path_buf())
    }

    /// The `total_source` of every mesh (photons/s).
    #[getter]
    fn total_source(&self) -> Option<f64> {
        self.file.total_source()
    }

    #[getter]
    fn ids(&self) -> Vec<u32> {
        self.file.tallies().map(|t| t.id).collect()
    }

    /// What each mesh is (no voxels read).
    #[getter]
    fn tallies(&self, py: Python<'_>) -> PyResult<Vec<PyCdgsInfo>> {
        self.file
            .tallies()
            .map(|t| PyCdgsInfo::new(py, t.clone()))
            .collect()
    }

    /// Meshes whose header could not be read, and why.
    #[getter]
    fn failed(&self) -> Vec<(u32, String)> {
        self.file.failed().to_vec()
    }

    fn info(&self, py: Python<'_>, id: u32) -> PyResult<PyCdgsInfo> {
        let info = self.file.tally(id).map_err(read_error)?.clone();
        PyCdgsInfo::new(py, info)
    }

    /// Mesh `id`'s source: per voxel and energy bin (or line), with a Total
    /// when there is more than one; `volume` per voxel.
    fn read(&self, py: Python<'_>, id: u32) -> PyResult<PyMeshTally> {
        let file = &self.file;
        let tally = py.allow_threads(|| file.read(id)).map_err(read_error)?;
        PyMeshTally::new(py, tally)
    }

    fn __getitem__(&self, py: Python<'_>, id: u32) -> PyResult<PyMeshTally> {
        self.read(py, id)
    }

    fn __len__(&self) -> usize {
        self.file.tallies().count()
    }

    fn __contains__(&self, id: u32) -> bool {
        self.file.tally(id).is_ok()
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<PyObject> {
        let ids = self.ids().into_pyobject(py)?;
        Ok(ids.call_method0("__iter__")?.unbind())
    }

    fn keys(&self) -> Vec<u32> {
        self.ids()
    }

    fn __repr__(&self) -> String {
        let path = self
            .file
            .path()
            .map_or("<memory>".to_string(), |p| p.display().to_string());
        format!("CdgsFile({path:?}, meshes={:?})", self.ids())
    }
}

/// What a mesh of a CDGS file is.
#[pyclass(name = "CdgsInfo", module = "aci_tallies.cdgs", frozen)]
pub struct PyCdgsInfo {
    info: CdgsInfo,
    #[pyo3(get)]
    geometry: Py<PyGeometry>,
    /// Energy bins (MeV), or gamma lines (`kind == "line"`).
    #[pyo3(get)]
    energy: Py<PyAxis>,
}

impl PyCdgsInfo {
    fn new(py: Python<'_>, info: CdgsInfo) -> PyResult<PyCdgsInfo> {
        Ok(PyCdgsInfo {
            geometry: geometry(py, &info.geometry)?,
            energy: Py::new(
                py,
                PyAxis {
                    axis: info.energy.clone(),
                },
            )?,
            info,
        })
    }
}

#[pymethods]
impl PyCdgsInfo {
    #[getter]
    fn id(&self) -> u32 {
        self.info.id
    }

    /// The line after `mesh_id`.
    #[getter]
    fn comment(&self) -> &str {
        &self.info.comment
    }

    /// Cooling time (s).
    #[getter]
    fn cooling_time(&self) -> f64 {
        self.info.cooling_time
    }

    /// The mesh's source strength (`mesh_source`, photons/s).
    #[getter]
    fn strength(&self) -> f64 {
        self.info.strength
    }

    #[getter]
    fn warnings(&self) -> Vec<PyReadWarning> {
        warnings(&self.info.warnings)
    }

    fn __repr__(&self) -> String {
        format!(
            "CdgsInfo(id={}, cooling_time={:e}, strength={:e})",
            self.info.id, self.info.cooling_time, self.info.strength
        )
    }
}
