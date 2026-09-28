//! Rust port of sampdist_1.f:
//! SAMPDIST - Sample distances in neighbourhoods
//! written by Mark Hill, January-June 2011
//!
//! For each location, writes the nearest `--neighbours` locations (itself
//! first) in increasing distance order. Distances are planar (easting,
//! northing; the original) or, with `--distance geodesic`, WGS84 geodesic
//! distances in km from (longitude, latitude).

use frescalo::cli::{fatal, Args, OptSpec, Session, NO_HOLD};
use frescalo::config::{limits, SampdistParams, DEFAULT_SPATIAL_NEIGHBOURS};
use frescalo::distance::{check_lon_lat, DistanceEngine, Metric};
use frescalo::fortran::{cout, name_to_string, DataReader};
use frescalo::neighbourhood::spatial::{looks_like_degrees, rank_by_distance, read_locations};
use frescalo::output::write_distance_row;
use std::io::{BufWriter, Write};

const PROGRAM: &str = "sampdist";

const OPTIONS: &[OptSpec] = &[
    OptSpec {
        name: "locations",
        value: Some("FILE"),
        help: "Locations file [site x y]",
    },
    OptSpec {
        name: "output",
        value: Some("FILE"),
        help: "Output file of neighbourhood distances (must not exist)",
    },
    OptSpec {
        name: "neighbours",
        value: Some("N"),
        help: "Nearest sites to keep per site, including itself [default: 200]",
    },
    OptSpec {
        name: "distance",
        value: Some("METHOD"),
        help: "planar (x y = easting northing) or geodesic (x y = lon lat, WGS84, km) [default: planar]",
    },
    NO_HOLD,
];

fn main() {
    let args = Args::parse(PROGRAM, "SAMPDIST - Sample distances in neighbourhoods", OPTIONS);
    let mut session = Session::new(PROGRAM, args.has("no-hold"));
    let metric = match args.get("distance") {
        None => Metric::Planar,
        Some(v) => Metric::from_option(v)
            .unwrap_or_else(|| fatal(PROGRAM, &format!("--distance must be 'planar' or 'geodesic', not '{}'", v))),
    };

    cout("");
    cout(" SAMPDIST - Sample distances in neighbourhoods");
    cout(" written by Mark Hill, January-June 2011");
    cout("");
    let (_filein, fin) = session.input_file(args.get("locations"), &[" Type name of file with locations ...."]);
    let mut reader = DataReader::new(fin);
    let (_fileou, fout) = session.output_file(
        args.get("output"),
        &[" Type name of output file with neighbourhood distances ..."],
    );
    let mut unit9 = BufWriter::new(fout);
    let params = SampdistParams {
        neighbours: session.neighbours(&args, DEFAULT_SPATIAL_NEIGHBOURS),
        metric,
    };

    let locs = read_locations(&mut reader, limits::SAMPDIST_SITES);

    match params.metric {
        Metric::Planar => {
            if args.get("distance").is_none() && locs.m > 1 && looks_like_degrees(&locs) {
                eprintln!(
                    "{}: warning: all coordinates look like longitude/latitude; planar distances on degrees \
                     distort neighbourhoods (consider --distance geodesic)",
                    PROGRAM
                );
            }
        }
        Metric::GeodesicWgs84 => {
            for i in 1..=locs.m {
                let c = &locs.coords[i];
                if let Err(e) = check_lon_lat(c.x64, c.y64) {
                    fatal(PROGRAM, &format!("site {}: {}", name_to_string(&locs.names[i]).trim(), e));
                }
            }
            cout(" Distances: geodesic on the WGS84 ellipsoid (Karney 2013), in km");
        }
    }
    if params.neighbours as usize > locs.m && args.get("neighbours").is_some() {
        eprintln!(
            "{}: warning: --neighbours {} exceeds the number of sites ({}); writing {} per site",
            PROGRAM, params.neighbours, locs.m, locs.m
        );
    }

    let engine = DistanceEngine::new(params.metric);
    let mut ranked = Vec::with_capacity(locs.m);
    let nout = (params.neighbours.max(0) as usize).min(locs.m);
    for i1 in 1..=locs.m {
        if i1 % 100 == 0 {
            frescalo::fortran::ld_line(&format!("Calculating distances   {}", frescalo::fortran::ld_i(locs.m as i64)));
        }
        rank_by_distance(i1, &locs, &engine, &mut ranked);
        for (is2, &(dist, iis2)) in ranked.iter().take(nout).enumerate() {
            write_distance_row(&mut unit9, &locs.names[i1], &locs.names[iis2 as usize], is2 as i64 + 1, dist as f32);
        }
    }

    unit9.flush().unwrap();
    session.finish();
}
