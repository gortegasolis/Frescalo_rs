//! Distances between sites: planar (the original) or geodesic.
//!
//! * **Planar** treats the two coordinates as easting and northing on a flat
//!   map (a national grid, UTM, …) and uses straight-line (Euclidean)
//!   distance, in the units of the input. This is the original behaviour and
//!   the default, computed in single precision exactly as the Fortran does.
//! * **Geodesic** treats the coordinates as longitude and latitude in degrees
//!   and measures the shortest distance over the WGS84 ellipsoid (Karney
//!   2013, via `geographiclib-rs`), in kilometres, in double precision. Use it
//!   when sites are given in degrees: Euclidean distance on degrees
//!   exaggerates east–west separation away from the equator (a degree of
//!   longitude shrinks with the cosine of latitude) and fails across the
//!   180° meridian and near the poles.
//!
//! Only the *rank* of a neighbour's distance is used downstream (by
//! `neighsim`'s spatial weight), so the choice of metric matters exactly when
//! it changes which sites are nearest, or their order.

use geographiclib_rs::{Geodesic, InverseGeodesic};

/// How distances between sites are measured.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Metric {
    /// Euclidean distance on projected easting/northing (input units).
    #[default]
    Planar,
    /// Ellipsoidal distance on WGS84 from longitude/latitude, in km.
    GeodesicWgs84,
}

impl Metric {
    /// Parse the `--distance` option value.
    pub fn from_option(s: &str) -> Option<Metric> {
        match s {
            "planar" => Some(Metric::Planar),
            "geodesic" => Some(Metric::GeodesicWgs84),
            _ => None,
        }
    }
}

/// Site coordinates: `(x, y)` = (easting, northing) in planar mode,
/// (longitude, latitude) in degrees in geodesic mode.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Coord {
    /// Planar-mode coordinates, parsed in single precision as in the original.
    pub x32: f32,
    pub y32: f32,
    /// Geodesic-mode coordinates, parsed in double precision.
    pub x64: f64,
    pub y64: f64,
}

/// Distance calculator, built once per run.
pub struct DistanceEngine {
    metric: Metric,
    geod: Option<Geodesic>,
}

impl DistanceEngine {
    pub fn new(metric: Metric) -> Self {
        DistanceEngine {
            metric,
            geod: match metric {
                Metric::Planar => None,
                Metric::GeodesicWgs84 => Some(Geodesic::wgs84()),
            },
        }
    }

    pub fn metric(&self) -> Metric {
        self.metric
    }

    /// Distance from `a` to `b`. Planar distances are computed in `f32` with
    /// the original operation order (then widened losslessly); geodesic
    /// distances are computed in `f64` and returned in kilometres.
    #[inline]
    pub fn distance(&self, a: &Coord, b: &Coord) -> f64 {
        match self.metric {
            Metric::Planar => {
                let de = a.x32 - b.x32;
                let dn = a.y32 - b.y32;
                (de * de + dn * dn).sqrt() as f64
            }
            Metric::GeodesicWgs84 => {
                let g = self.geod.as_ref().unwrap();
                //                      lat1   lon1   lat2   lon2
                let s12: f64 = g.inverse(a.y64, a.x64, b.y64, b.x64);
                s12 / 1000.0
            }
        }
    }
}

/// Check that a longitude/latitude pair is in range
/// (longitude −180…180, latitude −90…90).
pub fn check_lon_lat(lon: f64, lat: f64) -> Result<(), String> {
    if !(-180.0..=180.0).contains(&lon) {
        return Err(format!("longitude {} is outside -180 to 180", lon));
    }
    if !(-90.0..=90.0).contains(&lat) {
        return Err(format!("latitude {} is outside -90 to 90", lat));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn planar(x: f32, y: f32) -> Coord {
        Coord { x32: x, y32: y, ..Default::default() }
    }

    fn lonlat(lon: f64, lat: f64) -> Coord {
        Coord { x64: lon, y64: lat, ..Default::default() }
    }

    #[test]
    fn planar_345_triangle_zero_and_symmetry() {
        let e = DistanceEngine::new(Metric::Planar);
        assert_eq!(e.distance(&planar(0.0, 0.0), &planar(3.0, 4.0)), 5.0);
        assert_eq!(e.distance(&planar(7.0, 2.0), &planar(7.0, 2.0)), 0.0);
        let (a, b) = (planar(1.5, -2.0), planar(-4.0, 9.25));
        assert_eq!(e.distance(&a, &b), e.distance(&b, &a));
    }

    #[test]
    fn geodesic_known_distances() {
        let e = DistanceEngine::new(Metric::GeodesicWgs84);
        // One degree of longitude on the equator: a·π/180 with a = 6378.137 km.
        let d = e.distance(&lonlat(0.0, 0.0), &lonlat(1.0, 0.0));
        assert!((d - 111.319_490_793_273_57).abs() < 1e-9, "{}", d);
        // Pole to pole along a meridian: twice the WGS84 quarter meridian
        // (10 001.965 729 km).
        let d = e.distance(&lonlat(0.0, -90.0), &lonlat(0.0, 90.0));
        assert!((d - 20_003.931_458_6).abs() < 1e-6, "{}", d);
    }

    #[test]
    fn geodesic_argument_order() {
        // (lon 0, lat 60) to (lon 10, lat 60) is about 557 km; if longitude
        // and latitude were swapped it would be (lat 0 → lat 10 along lon 60),
        // about 1106 km.
        let e = DistanceEngine::new(Metric::GeodesicWgs84);
        let d = e.distance(&lonlat(0.0, 60.0), &lonlat(10.0, 60.0));
        assert!((d - 557.0).abs() < 2.0, "{}", d);
    }

    #[test]
    fn geodesic_crosses_antimeridian() {
        let e = DistanceEngine::new(Metric::GeodesicWgs84);
        let d = e.distance(&lonlat(179.5, 0.0), &lonlat(-179.5, 0.0));
        assert!((d - 111.319_490_793_273_57).abs() < 1e-6, "{}", d);
    }

    #[test]
    fn lon_lat_ranges() {
        assert!(check_lon_lat(-180.0, 90.0).is_ok());
        assert!(check_lon_lat(181.0, 0.0).is_err());
        assert!(check_lon_lat(0.0, -90.5).is_err());
    }
}
